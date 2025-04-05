%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SCPMinFuelSolver: An instance of SCPProblem class, for solving
% minimum-fuel trajectory optimization problems.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef SCPMinFuelSolver < SCPProblem
    properties
        p           % Class with all parameters (terminal conditions, etc.)
        init_guess_struct % Struct with initial guesses for each variable
        impose_trust_region_struct % Struct with boolean fields for imposing trust regions on optimization variables
        traj_with_stm % Trajectory with state transition matrix, obtained from the solution
        A
        B
        c
        d
        Nmaneuvers
    end

    methods
        function obj = SCPMinFuelSolver(p, init_guess_struct, impose_trust_region_struct)
            % Constructor for the SCPMinFuelSolver class
            arguments
                p % Class with all parameters (terminal conditions, etc.)
                init_guess_struct struct {has_x} % Struct with initial guesses for each variable, must have field x
                impose_trust_region_struct struct {has_x} = struct('x', true, 'u', false, 'v', false);
            end

            % Inherit from SCPProblem class
            obj@SCPProblem();

            obj.p = p;
            
            % Initialize the problem
            obj = obj.initialize(init_guess_struct, impose_trust_region_struct);
                
        end
        
        function obj = initialize(obj, init_guess_struct, impose_trust_region_struct)

            % Initialize the trust regions and initial guesses based on mesh_type and problem parameters           
            
            obj.init_guess_struct = init_guess_struct;
            obj.impose_trust_region_struct = impose_trust_region_struct;

            obj.trust_region_scaling = struct('x', 1, 'u', obj.p.u_max); % TODO: detailed if-else for trust region scaling

            if isfield(obj.init_guess_struct, 'v') && isfield(obj.impose_trust_region_struct, 'v') && obj.impose_trust_region_struct.v
                warning("Imposing a trust region on v (control magnitude) is not recommended. Consider using the trust region on u (control vector) instead");
            end 
            obj.impose_trust_region_struct.v = false;

            if obj.p.includeMass
                obj.p.nX = obj.p.DS.nx + 1; % orbital state + mass
                obj.p.nU = obj.p.DS.nu + 1; % control + control norm
                obj.init_guess_struct.x(obj.p.nX, :) = obj.p.Spacecraft.logmass_init * ones(obj.p.Nnode, 1); % mass variables
            end

            if isfield(obj.init_guess_struct, 'u')
                obj.init_guess_struct.v = vecnorm(obj.init_guess_struct.u);
            else
                switch obj.p.CntrlType
                    case 'LowThrust'
                        obj.init_guess_struct.u = zeros(obj.p.DS.nu, obj.p.Nseg);
                        obj.init_guess_struct.v = zeros(obj.p.Nseg, 1);
                    case 'impulsive'
                        obj.init_guess_struct.u = zeros(obj.p.DS.nu, obj.p.Nseg + 1);
                        obj.init_guess_struct.v = zeros(obj.p.Nseg + 1, 1);
                end
            end
            
            if strcmp(obj.p.mesh_type, 'adaptive')
                % Initialize according to the given time history as the initial guess
                if ~isfield(obj.init_guess_struct, 's')
                    obj.init_guess_struct.s = obj.p.Nseg * diff(obj.p.t_his);
                end
                obj.impose_trust_region_struct.s = 0;
                if strcmp(obj.p.CntrlType, 'LowThrust') && ~obj.p.includeMass
                    obj.init_guess_struct.J = v_ref0' * s_ref0 / obj.p.Nseg;
                    obj.impose_trust_region_struct.J = 1;
                end
            end

            if obj.p.is_parameterized_x_init
                obj.impose_trust_region_struct.tau_x_init = true;
                if ~isfield(obj.init_guess_struct, 'tau_x_init')
                    obj.init_guess_struct.tau_x_init = 0;
                end
            elseif isfield(obj.init_guess_struct, 'tau_x_init')
                warning('Initial guess has field tau_x_init, but is_parameterized_x_init is set to false');
            end
            if obj.p.is_parameterized_x_fin
                obj.impose_trust_region_struct.tau_x_fin = true;
                if ~isfield(obj.init_guess_struct, 'tau_x_fin')
                    obj.init_guess_struct.tau_x_fin = 0;
                end
            elseif isfield(obj.init_guess_struct, 'tau_x_fin')
                warning('Initial guess has field tau_x_fin, but is_parameterized_x_fin is set to false');
            end

            % Initialize the parent class
            initialize@SCPProblem(obj);

        end
            
        function vars = define_vars(obj)
            % DEFINE_VARS: Define the optimization variables
            vars.x = sdpvar(obj.p.DS.nx, obj.p.Nseg + 1, 'full');
            switch obj.p.CntrlType
                case 'LowThrust'
                    vars.u = sdpvar(obj.p.DS.nu, obj.p.Nseg, 'full');
                    vars.v = sdpvar(obj.p.Nseg, 1);
                case 'impulsive'
                    vars.u = sdpvar(obj.p.DS.nu, obj.p.Nseg + 1, 'full');
                    vars.v = sdpvar(obj.p.Nseg + 1, 1);
            end
            if strcmp(obj.p.mesh_type, 'adaptive')
                vars.s = sdpvar(obj.p.Nseg, 1);
            end
            if ~obj.p.includeMass && strcmp(obj.p.mesh_type, 'adaptive') && strcmp(obj.p.CntrlType, 'LowThrust')
                vars.J = sdpvar;
            end
            % Add variables for moving initial and final states
            if obj.p.is_parameterized_x_init
                vars.tau_x_init = sdpvar;
            end
            if obj.p.is_parameterized_x_fin
                vars.tau_x_fin = sdpvar;
            end
        end

        function J0 = objective(obj, vars)
            % OBJECTIVE: Evaluate the objective function
            % vars can be sdpvars or doubles, since the objective function is convex by assumption
            switch obj.p.CntrlType
                case 'impulsive'
                    J0 = sum(vars.v);
                case 'LowThrust'
                    if obj.p.includeMass
                        J0 = -vars.x(end, end);
                    else
                        switch obj.p.mesh_type
                            case 'fixed'
                                J0 = sum(vars.v);
                            case 'adaptive'
                                J0 = vars.J;
                        end
                    end
            end
        end
        
        function constraints = convex_eq(obj, vars)
            % CONVEX_EQ: Define the convex equality constraints
            % vars should be sdpvars. This funcition is called only once 
            constraints = [
                obj.initial_state_constraint(vars)
                obj.final_state_constraint(vars)
                obj.mass_flow_constraint(vars)
                obj.time_of_flight_constraint(vars)
            ];
        end

        function constraints = convex_ineq(obj, vars)
            % CONVEX_INEQ: Define the convex inequality constraints
            % vars should be sdpvars. This funcition is called only once
            constraints = [
                obj.control_cone_constraint(vars)
                obj.control_bounds_constraint(vars)
                obj.logmass_constraint(vars)
                obj.decreasing_logmass_constraint(vars)
                obj.s_bounds_constraint(vars)
                obj.J_constraint(vars)
            ];
        end

        function constraintLHS = noncvx_eq(obj, vars)
            % NONCVX_EQ: Define LHS of nonconvex equality constraints
            % vars should be doubles
            constraintLHS = [
                obj.nonlin_eom_constraint(vars)
                obj.moving_initial_state_constraint(vars)
                obj.moving_final_state_constraint(vars)
            ];

        end
        
        function constraintLHS = noncvx_ineq(obj, vars)
            % NONCVX_INEQ: Define LHS of nonconvex inequality constraints
            constraintLHS = [];
        end

        function constraintLHS = noncvx_eq_relaxed(obj, vars, ref_vars)
            constraintLHS = [
                obj.linearized_eom(vars, ref_vars)
                obj.moving_initial_state_constraint_linearized(vars, ref_vars)
                obj.moving_final_state_constraint_linearized(vars, ref_vars)
            ];
                
        end

        function constraintLHS = noncvx_ineq_relaxed(obj, vars, ref_vars)
            constraintLHS = [];
        end
    
        function constraint = initial_state_constraint(obj, vars)
            % Define the initial state constraint
            if obj.p.is_parameterized_x_init
                constraint = [];
            else
                switch obj.p.CntrlType
                    case 'LowThrust'
                        constraint = vars.x(1:obj.p.DS.nx, 1) - obj.p.x_init == 0;
                    case 'impulsive'
                        constraint = vars.x(1:obj.p.DS.nx, 1) - obj.p.DS.add_impulse(obj.p.x_init, vars.u(1:obj.p.DS.nu, 1)) == 0;
                end
                constraint = [constraint:'initial state'];
            end
        end

        function constraint = final_state_constraint(obj, vars)
            % Define the final state constraint
            if obj.p.is_parameterized_x_fin
                constraint = [];
            else            
                constraint = [vars.x(1:obj.p.DS.nx, end) - obj.p.x_fin == 0]: 'final state';
            end
        end

        function constraint = mass_flow_constraint(obj, vars)
            % Define the mass flow constraint
            if obj.p.includeMass
                logmass = vars.x(obj.p.nX, :);
                Gamma = vars.u(obj.p.nU, :)';
                constraint = [];
                for k = 1:obj.p.Nseg
                    constraint = [constraint; logmass(k + 1) - logmass(k) + ...
                        (obj.p.t_his(k + 1) - obj.p.t_his(k)) / obj.p.ex_vel * Gamma(k) == 0];
                end
                constraint = [constraint: 'mass flow'];
            else
                constraint = [];
            end
        end

        function constraint = time_of_flight_constraint(obj, vars)
            % Define the time-of-flight constraint for adaptive mesh
            if strcmp(obj.p.mesh_type, 'adaptive')
                constraint = [sum(vars.s) / obj.p.Nseg - obj.p.tof == 0]: 'TOF';
            else
                constraint = [];
            end
        end   

        function constraint = control_cone_constraint(obj, vars)
            % Define the control cone constraint
            if obj.p.includeMass
                Gamma = vars.u(obj.p.nU, :)';
                constraint = cone([Gamma'; vars.u(1:obj.p.DS.nu, :)]);
            else
                constraint = cone([vars.v'; vars.u(1:obj.p.DS.nu, :)]);
            end
            constraint = [constraint: 'control cone'];
        end

        function constraint = control_bounds_constraint(obj, vars)
            % Define the lower/upper bounds on the control magnitude
            if obj.p.includeMass
                error('Control bounds constraint not yet implemented for mass flow');
            elseif strcmp(obj.p.CntrlType, 'LowThrust')
                constraint = [vars.v <= obj.p.u_max];
                % constraint = [vars.v <= 1];
                if ~isempty(obj.p.u_min)
                    constraint = [constraint; obj.p.u_min <= vars.v];
                    % constraint = [constraint; obj.p.u_min / obj.p.u_max <= vars.v];
                end
                constraint = [constraint: 'control bounds'];
            else
                constraint = [];
            end
        end

        function constraint = logmass_constraint(obj, vars)
            % Define the log-mass constraint
            if obj.p.includeMass
                constraint = [Gamma - (obj.p.T_max * exp(-obj.ref_vars.x(obj.p.nX, :))) <= 0]: 'log-mass';
            else
                constraint = [];
            end
        end

        function constraint = decreasing_logmass_constraint(obj, vars)
            % Define the decreasing log-mass constraint
            if obj.p.includeMass
                logmass = vars.x(obj.p.nX, :);
                constraint = [logmass(:) - obj.p.logmass_init <= 0]: 'decreasing log-mass';
            else
                constraint = [];
            end
        end

        function constraint = s_bounds_constraint(obj, vars)
            % Define the bounds on the adaptive mesh parameter s
            if strcmp(obj.p.mesh_type, 'adaptive')
                constraint = obj.p.s_min - vars.s <= 0;
                if ~isempty(obj.p.s_max)
                    constraint = [constraint; vars.s - obj.p.s_max <= 0];
                end
                constraint = [constraint: 's bounds'];
            else
                constraint = [];
            end
        end

        function constraint = J_constraint(obj, vars)
            % Define the constraint on the objective function J for adaptive mesh
            if ~obj.p.includeMass && strcmp(obj.p.mesh_type, 'adaptive') && strcmp(obj.p.CntrlType, 'LowThrust')
                constraint = [((vars.v - obj.ref_vars.v)' * obj.ref_vars.s ...
                    + obj.ref_vars.v' * (vars.s - obj.ref_vars.s)) / obj.p.Nseg - vars.J <= 0]: 'J constraint';
            else
                constraint = [];
            end
        end

        function constraintLHS = moving_initial_state_constraint(obj, vars)
            if obj.p.is_parameterized_x_init
                switch obj.p.CntrlType
                    case 'LowThrust'
                        constraintLHS = vars.x(1:obj.p.DS.nx, 1) - obj.p.eval_x_init(vars.tau_x_init);
                    case 'impulsive'
                        constraintLHS = vars.x(1:obj.p.DS.nx, 1) ...
                        - obj.p.DS.add_impulse(obj.p.eval_x_init(vars.tau_x_init), vars.u(1:obj.p.DS.nu, 1));
                end
            else
                constraintLHS = [];
            end
        end

        function constraintLHS = moving_initial_state_constraint_linearized(obj, vars, ref_vars)
            if obj.p.is_parameterized_x_init
                switch obj.p.CntrlType
                    case 'LowThrust'
                        constraintLHS = vars.x(1:obj.p.DS.nx, 1) ...
                                    - obj.p.eval_x_init(ref_vars.tau_x_init)...
                                    - obj.p.eval_x_init_gradient(ref_vars.tau_x_init) * (vars.tau_x_init - ref_vars.tau_x_init);
                    case 'impulsive'
                        constraintLHS = vars.x(1:obj.p.DS.nx, 1) ...
                                    - obj.p.DS.add_impulse(obj.p.eval_x_init(ref_vars.tau_x_init), vars.u(1:obj.p.DS.nu, 1))...
                                    - obj.p.eval_x_init_gradient(ref_vars.tau_x_init) * (vars.tau_x_init - ref_vars.tau_x_init);
                end
            else
                constraintLHS = [];
            end
        end

        function constraintLHS = moving_final_state_constraint(obj, vars)
            if obj.p.is_parameterized_x_fin
                constraintLHS = vars.x(1:obj.p.DS.nx, end) - obj.p.eval_x_fin(vars.tau_x_fin);
            else
                constraintLHS = [];
            end
        end

        function constraintLHS = moving_final_state_constraint_linearized(obj, vars, ref_vars)
            if obj.p.is_parameterized_x_fin
                constraintLHS = vars.x(1:obj.p.DS.nx, end) ...
                            - obj.p.eval_x_fin(ref_vars.tau_x_fin)...
                            - obj.p.eval_x_fin_gradient(ref_vars.tau_x_fin) * (vars.tau_x_fin - ref_vars.tau_x_fin);
            else
                constraintLHS = [];
            end
        end

         function constraintLHS = obstacle_avoidance_constraint(obj, vars)
            % OBSTACLE_AVOIDANCE_CONSTRAINT: Define the constraint for avoiding obstacles
            constraintLHS = [];
            if obj.p.n_obs == 0
                return;
            end
            for k = 1:obj.p.Nseg
                p_k = vars.x(1:3, k);
                for i = 1:obj.p.n_obs
                    p_obs = obj.p.p_obs(:, i);
                    R_obs = obj.p.R_obs(i);
                    constraintLHS = [constraintLHS
                        R_obs - norm(p_k - p_obs)];
                end
            end
        end

        function constraintLHS = obstacle_avoidance_constraint_linearized(obj, vars, ref_vars)
            % OBSTACLE_AVOIDANCE_CONSTRAINT_LINEARIZED: Linearize the constraint for avoiding obstacles
            constraintLHS = [];
            if obj.p.n_obs == 0
                return;
            end
            for k = 1:obj.p.Nseg
                p_k = vars.x(1:3, k);
                p_k_ref = ref_vars.x(1:3, k);
                for i = 1:obj.p.n_obs
                    p_obs = obj.p.p_obs(:, i);
                    R_obs = obj.p.R_obs(i);
                    constraintLHS = [constraintLHS
                        R_obs - norm(p_k_ref - p_obs) - unit(p_k_ref - p_obs)' * (p_k - p_k_ref)];
                end
            end
        end

        function constraintLHS = nonlin_eom_constraint(obj, vars)
            constraintLHS = [];
            t_his = obj.get_t_his(vars);
            switch obj.p.CntrlType
                case 'LowThrust'
                    for k = 1:obj.p.Nseg
                        [~, x] = obj.p.DS.propagate_with_LT(vars.x(1:obj.p.DS.nx, k), t_his(k:k+1), vars.u(1:obj.p.DS.nu, k));
                        x_kplus1 = x(end, 1:obj.p.DS.nx)';
                        constraintLHS = [constraintLHS; vars.x(1:obj.p.DS.nx, k+1) - x_kplus1];
                    end
                case 'impulsive'
                    for k = 1:obj.p.Nseg
                        [~, x] = obj.p.DS.propagate(vars.x(1:obj.p.DS.nx, k), t_his(k:k+1));
                        x_kplus1_preDV = x(end, 1:obj.p.DS.nx)';
                        x_kplus1 = obj.p.DS.add_impulse(x_kplus1_preDV, vars.u(1:obj.p.DS.nu, k+1));
                        constraintLHS = [constraintLHS; vars.x(1:obj.p.DS.nx, k+1) - x_kplus1];
                    end
            end
        end
        

        function constraintLHS = linearized_eom(obj, vars, ref_vars)
        %% LINEARIZED_EOM: Construct the LHS of the linear time-varying dynamics constraint

            nx = obj.p.DS.nx;
            nu = obj.p.DS.nu;
            Nseg = obj.p.Nseg;

            constraintLHS = [];

            if strcmp(obj.p.CntrlType, 'LowThrust') && strcmp(obj.p.mesh_type, 'fixed')
                [obj.A, obj.B, obj.c] = obj.p.DS.discretize_LT(ref_vars.x(1:nx,:), ref_vars.u(1:nu,:), obj.p.t_his);
                for k = 1:Nseg
                    constraintLHS = [constraintLHS; vars.x(1:nx,k+1) - obj.A(1:nx,1:nx,k) * vars.x(1:nx,k) -...
                                     obj.B(1:nx, 1:nu,k)* vars.u(1:nu,k)- obj.c(1:nx,:,k)];
                end
            elseif strcmp(obj.p.CntrlType, 'LowThrust') && strcmp(obj.p.mesh_type, 'adaptive')
                tau_his = linspace(0, 1, Nseg+1);
                [obj.A, obj.B, obj.c, obj.d] = obj.p.DS.discretize_LT_adaptive(ref_vars.x(1:nx,:), ref_vars.u(1:nu,:), ref_vars.s, tau_his);
                for k = 1:Nseg
                    constraintLHS = [constraintLHS; vars.x(1:nx,k+1) - obj.A(1:nx,1:nx,k) * vars.x(1:nx,k) - ...
                                        obj.B(1:nx, 1:nu,k)*vars.u(1:nu,k) - obj.c(1:nx,:,k) - obj.d(:,:,k)* vars.s(k)];
                end
            elseif strcmp(obj.p.CntrlType, 'impulsive') && strcmp(obj.p.mesh_type, 'fixed')
                [obj.A, obj.B, obj.c] = obj.p.DS.discretize_impulsive(ref_vars.x(1:nx,:), obj.p.t_his);
                for k = 1:Nseg
                    constraintLHS = [constraintLHS; vars.x(1:nx,k+1) - obj.A(1:nx,1:nx,k) * vars.x(1:nx,k) - ...
                                         obj.B(1:nx, 1:nu,k+1)*vars.u(1:nu,k+1) - obj.c(1:nx,:,k)];
                end % Note the indexing of vars.u
            elseif strcmp(obj.p.CntrlType, 'impulsive') && strcmp(obj.p.mesh_type, 'adaptive')
                tau_his = linspace(0, 1, Nseg+1);
                [obj.A, obj.B, obj.c, obj.d] = obj.p.DS.discretize_impulsive_adaptive(ref_vars.x(1:nx,:), ref_vars.s, tau_his);
                for k = 1:Nseg
                    constraintLHS = [constraintLHS; vars.x(1:nx,k+1) - obj.A(1:nx,1:nx,k) * vars.x(1:nx,k) - ...
                                     obj.B(1:nx, 1:nu,k+1)*vars.u(1:nu,k+1) - obj.c(1:nx,:,k) - obj.d(:,:,k)* vars.s(k)];
                end % Note the indexing of vars.u
            else
                error('Invalid control type or mesh type');
            end
        end
            
        function t_his = get_t_his(obj, sol)
            % GET_T_HIS: Get the time history from the solution
            % sol should be doubles from the solution of the optimization problem
            arguments
                obj
                sol = obj.sol
            end

            switch obj.p.mesh_type
                case 'fixed'
                    t_his = obj.p.t_his;
                case 'adaptive'
                    t_his = [0; cumsum(sol.s)/obj.p.Nseg];
            end
        end
        

        %% Post Processing functions
        function set_traj_with_stm(obj, sol)
            % SET_TRAJ_WITH_STM: Get the trajectory with the state transition matrix from the solution
            % sol should be doubles from the solution of the optimization problem
            arguments
                obj
                sol = obj.sol
            end
            p = obj.p;

            t_his = obj.get_t_his(sol);

            if p.is_parameterized_x_init
                x_init = p.eval_x_init(sol.tau_x_init);
            else
                x_init = p.x_init;
            end
            
            switch obj.p.CntrlType
                case 'impulsive'
                    x_kplus1 = p.DS.add_impulse(x_init, sol.u(1:p.DS.nu, 1));
                    traj_with_stm = p.DS.propagate_with_STM(x_kplus1, t_his(1:2));

                    for k = 2:p.Nseg
                        x_kplus1_preDV = traj_with_stm.y(:, end);
                        x_kplus1 = p.DS.add_impulse(x_kplus1_preDV, sol.u(1:p.DS.nu, k));
                        traj_with_stm = odextend(traj_with_stm, [], t_his(k+1), x_kplus1);
                    end
                case 'LowThrust'
                    % u_handle = @(t) getZOH(t, t_his, sol.u);
                    % traj_with_stm = p.DS.propagate_with_LT_with_STM(x_init, t_his, u_handle);

                    % If the convergence tolerance is not small enough, small infeasibilities in the
                    % dynamics will result in not reach the target. 
                    traj_with_stm = p.DS.propagate_with_LT_with_STM(x_init, t_his(1:2), sol.u(1:p.DS.nu, 1));
                    for k = 2:p.Nseg
                        traj_with_stm = odextend(traj_with_stm, ...
                                                @(t,y) p.DS.EOM_with_LT_with_STM(t,y,sol.u(1:p.DS.nu, k)), ...
                                                t_his(k+1),...
                                                [obj.sol.x(1:p.DS.nx, k); vec(traj_with_stm.y(p.DS.nx+1:end, end))]);
                    end
            end
            obj.traj_with_stm = traj_with_stm;
        end

        function plot_trajectory(obj, options)
            arguments
                obj
                options.plot_maneuvers = true
            end
            figure; view(3)
            if isempty(obj.traj_with_stm)
                obj.set_traj_with_stm();
            end
            hold on 
            plot3d(obj.traj_with_stm.y(1:3,:), 'k', 'HandleVisibility', 'off');
            if options.plot_maneuvers
                quiver3d(obj.sol.x(:,1:size(obj.sol.u,2)), obj.sol.u, 'r', 'HandleVisibility', 'off');
            end
        end

        function plot_control_magnitude(obj, options)
            arguments
                obj
                options.dimensionalize = true
            end

            t_his = obj.get_t_his(obj.sol);
            u_norms = obj.sol.v;
            if options.dimensionalize
                t_his = t_his * obj.p.DS.SF.TU2DAYS;
                u_norms = u_norms * obj.p.DS.SF.a * 1E6;
            end
            figure
            switch obj.p.CntrlType
                case 'LowThrust'
                    stairsZOH(t_his, u_norms, 'k');
                case 'impulsive'
                    warning('Need to implement dimensionalization for impulsive control type');
                    stem(t_his, obj.sol.v, 'k', 'filled', 'LineWidth', 1.5);
            end
            
            if options.dimensionalize
                xlabel('$t$ [days]')
                ylabel('$\|u\|$ [mm/s$^2$]')
            else
                xlabel('$t$ [TU]')
                ylabel('$\|u\|$')
            end
        end

        function varargout = plot_primer_vector(obj)
        %% PLOT_PRIMER_VECTOR: Plot the primer vector magnitude from the solution
            if strcmp(obj.p.CntrlType, 'LowThrust')
                error('Primer vector plot is not implemented for LowThrust control type');
            end
            if isempty(obj.traj_with_stm)
                obj.set_traj_with_stm();
            end
            % Get the primer vector from the trajectory
            t_his = obj.get_t_his(obj.sol);
            t_eval = [];
            p_eval = [];
            for k = 1:obj.p.Nseg
                t1 = t_his(k);
                t2 = t_his(k+1);
                DV1 = obj.sol.u(:,k);
                DV2 = obj.sol.u(:,k+1);
                
                [t_eval_k, p_eval_k] = get_primer_vector_from_traj(obj.traj_with_stm, t1, t2, DV1, DV2);
                t_eval = [t_eval, t_eval_k];
                p_eval = [p_eval, p_eval_k];
            end
            pmag_eval = vecnorm(p_eval);
            % Identify peaks in the primer vector magnitudes
            [is_peak] = islocalmax(pmag_eval);
            % filter out the peaks that have primer magnitude close to unity
            is_peak = is_peak & pmag_eval > 1.0001;
            pmag_peaks = pmag_eval(is_peak);
            t_peaks = t_eval(is_peak);
            % Plot the primer vector magnitude
            figure
            hold on
            plot(t_eval, pmag_eval, 'k', HandleVisibility='off')
            plot(t_his, ones(size(t_his)), 'bo', DisplayName="$\Delta V$'s")
            plot(t_peaks, pmag_peaks, 'r*', DisplayName='Peaks')
            xlabel('$t$ [TU]')
            ylabel('$\|p\|$')
            legend()

            if nargout > 0
                varargout{1} = t_peaks;
                varargout{2} = pmag_peaks;
            end
        end

        function print_solution(obj)
            % PRINT_SOLUTION: Print the solution information
            assert(strcmp(obj.p.CntrlType, 'impulsive'), 'print_solution not implemented for LowThrust control type');
            disp(' ')
            disp('-----Solution Information-----')
            disp("Total Delta V: " + sum(obj.sol.v) * obj.p.DS.SF.v * 1000 + " m/s")
            t_his_days = obj.get_t_his(obj.sol) * obj.p.DS.SF.TU2DAYS;
            for k = 1:obj.p.Nseg + 1
                disp("DV" + k + ": " + obj.sol.v(k) * obj.p.DS.SF.v * 1000 + " m/s; " + ...
                    t_his_days(k) + " days")
            end
            disp("Total TOF: " + t_his_days(end) + " days")
            disp('-------------------------------')
            disp(' ')
        end                    

    end

end

%% Helper functions
function has_x(init_guess_struct)
    assert(isfield(init_guess_struct, 'x'), 'Initial guess must have the field x')
end

function [t_eval, p_eval] = get_primer_vector_from_traj(traj_with_stm, t1, t2, DV1, DV2, t_eval)
    %% get_primer_vector_from_traj: get the continuous primer vector history
    % from two reference burns
    % traj_with_stm: struct with augmented state with initial condition as identity matrix
    % at time t_0, defined from time t_0 to t_f
    % t1: time of first burn
    % t2: time of second burn
    % DV1: delta-v vector at t1
    % DV2: delta-v vector at t2
    % t_eval: time vector to evaluate the primer vector
    % p_eval: primer vector at times t_eval

    if nargin < 6
        t_eval = linspace(t1, t2, 100);
    end
    p1 = DV1/norm(DV1);
    p2 = DV2/norm(DV2);

    num_eval_pts = length(t_eval);
    stm_t_t0 = deval(traj_with_stm, t_eval, 7:42);
    stm_t_t0 = reshape(stm_t_t0, 6, 6, []);
    stm_t1_t0 = reshape(deval(traj_with_stm, t1, 7:42), [6 6]);
    stm_t2_t0 = reshape(deval(traj_with_stm, t2, 7:42), [6 6]);
    stm_t_t1 = pagemrdivide(stm_t_t0, stm_t1_t0);
    stm_t2_t1 = stm_t2_t0 / stm_t1_t0;
    M_t_t1 = stm_t_t1(1:3,1:3,:);
    N_t_t1 = stm_t_t1(1:3,4:6,:);
    M_t2_t1 = stm_t2_t1(1:3,1:3);
    N_t2_t1 = stm_t2_t1(1:3,4:6);
    N_t_t2 = pagemrdivide(N_t_t1, N_t2_t1);
    p_eval = NaN(length(p1), num_eval_pts);
    for i = 1:num_eval_pts
        p_eval(:,i) = (M_t_t1(:,:,i) - N_t_t2(:,:,i)*M_t2_t1)*p1 ...
            + N_t_t2(:,:,i) * p2;
    end
end
