classdef ExampleClass2 < SCPProblem
    
    properties
        init_guess_struct
        impose_trust_region_struct = struct('x', true, 'u', true, 'Gamma', false);
    end

    properties (Constant)
        g = [-9.81; 0; 0]; % m/s^2, gravity
        k_D = 0.5; % drag coefficient
        R_obs = [1; 1]; % m, radius of obstacles
        p_obs = [0 3 0.45; 0 7 -0.45]'; % m, position of obstacles
        n_obs = 2; % number of obstacles
        t_N = 5; % s
        x_init = [0 0 0 0 0.5 0]'; % m, m/s, initial state
        x_fin = [0 10 0 0 0.5 0]'; % m, m/s, final state
        m = 0.3; % kg, mass of quadrotor
        T_min = 1.0; % N, minimum thrust
        T_max = 4.0; % N, maximum thrust
        theta_max = pi/4; % rad, maximum tilt angle
        Nseg = 30; % number of segments
    end

    properties
        DS
        DeltaT
        t_his
    end

    properties
        A
        B
        c
    end

    methods

        function obj = ExampleClass2()
            % Must inherit from SCPProblem class
            obj@SCPProblem();

            obj.DS     = Quadrotor(obj.m, obj.k_D, obj.g);
            obj.DeltaT = obj.t_N/(obj.Nseg);
            obj.t_his  = linspace(0, obj.t_N, obj.Nseg+1);

            % Set the initial guess
            obj.init_guess_struct.x     = linspace_vec(obj.x_init, obj.x_fin, obj.Nseg+1);
            obj.init_guess_struct.u     = -obj.m * obj.g .* ones(obj.DS.nu, obj.Nseg);
            obj.init_guess_struct.Gamma = norm(obj.m * obj.g) * ones(1,obj.Nseg);

            % Must initialize the parent class
            obj.initialize()
        end

        function vars = define_vars(obj)
            vars.x = sdpvar(obj.DS.nx, obj.Nseg+1);
            vars.u = sdpvar(obj.DS.nu, obj.Nseg);
            vars.Gamma = sdpvar(1, obj.Nseg);
        end

        function J0 = objective(obj, vars)
            J0 = sum(vars.Gamma) * obj.DeltaT;
        end

        function constraints = convex_eq(obj, vars)
            constraints = [
                vars.x(:,1)   == obj.x_init
                vars.x(:,end) == obj.x_fin
                vars.u(:,1)   == -obj.m * obj.g
                vars.u(:,end) == -obj.m * obj.g
                vars.x(1,:)   == 0
            ];
        end

        function constraints = convex_ineq(obj, vars)
            constraints = [
                cone([vars.Gamma; vars.u])
                obj.T_min <= vars.Gamma
                vars.Gamma <= obj.T_max
                cos(obj.theta_max) * vars.Gamma <= vars.u(1,:)
            ];
        end

        function update_parameters(obj, ref_vars)
            [obj.A, obj.B, obj.c] = obj.DS.discretize_LT(ref_vars.x, ref_vars.u, obj.t_his);
        end

        function constraintLHS = noncvx_eq(obj, vars)
            constraintLHS = [];
            % Dynamics equation
            for k = 1:obj.Nseg
                [~, x] = obj.DS.propagate_with_LT(vars.x(:,k), obj.t_his(k:k+1), vars.u(:,k));
                x_kplus1 = x(end,:)';
                constraintLHS = [constraintLHS; vars.x(:,k+1) - x_kplus1];
            end
        end

        function constraintLHS = noncvx_ineq(obj, vars)
            constraintLHS = [];
            % Collision avoidance constraint
            for k = 1:obj.Nseg+1
                for i = 1:obj.n_obs
                    constraintLHS = [constraintLHS
                        obj.R_obs(i) - norm(vars.x(1:3,k) - obj.p_obs(:,i))];
                end
            end
        end

        function constraintLHS = noncvx_eq_relaxed(obj, vars, ref_vars)
            constraintLHS = [];

            for k = 1:obj.Nseg
                constraintLHS = [constraintLHS
                vars.x(:,k+1) - obj.A(:,:,k) * vars.x(:,k) - obj.B(:,:,k)* vars.u(:,k) - obj.c(:,:,k)];
            end
        end

        function constraintLHS = noncvx_ineq_relaxed(obj, vars, ref_vars)
            constraintLHS = [];
            for k = 1:obj.Nseg+1
                for i = 1:obj.n_obs
                    % Linearize the constraint
                    constraintLHS = [constraintLHS
                        obj.R_obs(i) - norm(ref_vars.x(1:3,k) - obj.p_obs(:,i)) ...
                        - unit(ref_vars.x(1:3,k) - obj.p_obs(:,i))' * (vars.x(1:3,k) - ref_vars.x(1:3,k))
                    ];
                end
            end
        end

    end
end
