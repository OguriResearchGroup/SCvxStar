%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SCvxStar: A SCP solver with trust region and augmented Lagrangian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef SCvxStar < handle
    properties
        this_iter SCPIter
        prev_iter SCPIter
        report = struct('iters', [], 'time', [], 'converged', [])
        iter_hist = SCPIter.empty % History of iterations
        clear_every_iter % Clear YALMIP variables every iteration
    end

    properties (SetAccess = immutable)
        scp_prob      % Instance of the SCPProblem class
        starttime     % Timer for optimization process
        constParams   % Structure with constant parameters, see SCPParams.m
        use_optimizer % Boolean to use YALMIP optimizer
        optimizer     % YALMIP optimizer object
        var_names     % names of the optimization variables
    end
    
    methods
        % Constructor
        function obj = SCvxStar(scp_prob, constParams, use_optimizer)
            arguments
                scp_prob SCPProblem
                constParams = SCPParams()
                use_optimizer = false
            end

            obj.starttime          = datetime('now');
            obj.scp_prob           = scp_prob;
            obj.use_optimizer      = use_optimizer;
            obj.constParams        = constParams;
            obj.this_iter          = SCPIter(1, obj);
            obj.this_iter.ref_vars = scp_prob.init_guess_struct;
            obj.this_iter.g_ref    = scp_prob.noncvx_eq(obj.this_iter.ref_vars);
            obj.this_iter.h_ref    = scp_prob.noncvx_ineq(obj.this_iter.ref_vars);
            obj.this_iter.J0_ref   = scp_prob.objective(obj.this_iter.ref_vars);
            obj.var_names          = fieldnames(scp_prob.sdp_vars);

            if obj.use_optimizer
                obj.optimizer = obj.create_optimizer();
            end
        end
        
        % Main function to solve the optimization problem
        function [converged] = solve(obj, options)
            arguments
                obj
                options.verbose = true
                options.clear_every_iter = false % Clear YALMIP variables every iteration
                options.resume = false % Resume from a specific iteration
            end % TODO: add levels of storing iterations: all, only necessary

            verbose = options.verbose;

            if options.resume
                iter = input('Resume from iteration: ');
                assert(iter > 1, 'Resume iteration must be greater than 1');
                assert(iter <= length(obj.iter_hist), 'Resume iteration must be less than the number of iterations already run');
                obj.prev_iter = obj.iter_hist(iter-1);
                obj.this_iter = obj.iter_hist(iter);
            else
                iter = 1;
                
            end

            converged = false;
            last_checked_close_to_optimum = 0;

            while iter < obj.constParams.k_max

                obj.scp_prob.pre_iteration();

                if obj.use_optimizer
                    obj.this_iter.solve_flag = obj.solve_subproblem_optimizer();
                else
                    obj.this_iter.solve_flag = obj.solve_subproblem();
                end

                % infeasible or other issues
                if obj.this_iter.solve_flag ~=0 && obj.this_iter.solve_flag ~=4 
                    keyboard
                    obj.create_report(iter, false);
                    return;
                end

                % Post iteration routine defined in the problem class
                obj.scp_prob.post_iteration();

                tic;
                [next_iter] = obj.postprocess(verbose, options.clear_every_iter);
                obj.this_iter.time_postprocessing = toc;

                obj.iter_hist(end+1) = obj.this_iter;

                obj.print_iteration_info(verbose);

                if obj.this_iter.solve_flag == 0
            
                    converged = obj.check_stopping_criteria();

                    if converged
                        obj.create_report(iter, true);
                        obj.fprintf_verbose(verbose, '\nConverged in %i iterations in %s\n', iter, obj.report.time);
                        return;
                    end

                    stalled = obj.check_stall(iter);

                    if stalled
                        keyboard
                        if obj.this_iter.chi > 1E-3
                            obj.fprintf_verbose(verbose, '\nAlgorithm is stalled. Exiting...\n');
                            obj.create_report(iter, false);
                            return;
                        else
                            fprintf('\nStalled but close to convergence. The trust region might have shrunk rapidly due to incorrect constraint linearization, etc.\n');
                            % Decrease penalty and enlarge trust region.)
                            % next_iter.r = next_iter.r * 10;
                            % next_iter.w = next_iter.w / 10;
                        end
                    end

                end

                if obj.check_close_to_optimum(iter, last_checked_close_to_optimum)
                    % beep;
                    disp('Recent iterations satisfy feasibility tolerance but objective function seems to be jumping around.')
                    % user_input = input('\nTo accept the best solution found so far, input y. Press any other key to continue: ', 's');
                    % switch user_input
                    %     case 'y'
                    %         % find the iteration with feasible solution and minimum objective function
                    %         feasible_iters = find([obj.iter_hist.chi] < obj.constParams.tol_feas);
                    %         [~, min_obj_idx] = min([obj.iter_hist(feasible_iters).J0]);
                    %         min_obj_iter = feasible_iters(min_obj_idx);                       
                    %         obj.create_report(min_obj_iter, true);
                    %         obj.fprintf_verbose(verbose, '\nAccepting solution from iteration %i and exiting.\n', min_obj_iter);
                    %         return;
                    %     otherwise
                    %         last_checked_close_to_optimum = iter;
                    % end
                end

                if obj.check_numerical(iter)
                    obj.fprintf_verbose(verbose, '\nNumerical issues detected. Exiting...\n');
                    obj.create_report(iter, false);
                    return;
                end          

                iter = iter + 1;
                obj.prev_iter = obj.this_iter;
                obj.this_iter = next_iter;
            end

            obj.fprintf_verbose(verbose, '\nReached %i iterations. To continue, input the number of additional iterations.', obj.constParams.k_max);
            beep;

            % Get user input for the number of additional iterations
            additional_iters = input('Additional iterations: ');
            if isempty(additional_iters)
                obj.create_report(iter, false);
                return;
            end

            while ~isnumeric(additional_iters) || additional_iters < 0
                disp('Input a valid number');
                additional_iters = input('Additional iterations: ');
            end

            obj.constParams.k_max = obj.constParams.k_max + additional_iters;

        end          

        % Solves the subproblem for the current iteration
        function [solve_flag] = solve_subproblem(obj)
            % Get the sdp object of the difference from the reference

            DeltaZ = obj.getDeltaZ(obj.scp_prob.sdp_vars);
            trust_region_constraint_lhs =  obj.get_trust_region_constraint_lhs(obj.scp_prob.sdp_vars);

            % Stack all constraints
            constraints_all = [
                % affine equality constraints
                obj.scp_prob.sdp_convex_eq
                % convex inequality constraints
                obj.scp_prob.sdp_convex_ineq
                % convexified nonconvex equality constraints
                obj.scp_prob.sdp_noncvx_eq_relaxed
                % convexified nonconvex inequality constraints
                obj.scp_prob.sdp_noncvx_ineq_relaxed
                % convexified nonconvex equality constraints that are imposed without slack variables
                obj.scp_prob.sdp_convexified_exact
                % nonconvex inequality slack variables are positive
                obj.scp_prob.slack_noncvx_ineq >= 0
                % trust region constraint
                trust_region_constraint_lhs <= 0
            ];

            J0 = obj.scp_prob.sdp_objective;
            P_cvx = obj.penalty_function(obj.scp_prob.slack_noncvx_eq, ...
                                        obj.scp_prob.slack_noncvx_ineq, ...
                                        DeltaZ=DeltaZ, ...
                                        trust_region_constraint_lhs=trust_region_constraint_lhs);

            % Augment the objective function with the penalty term
            L = J0 + P_cvx;

            % Solve the subproblem
            [diagnostics, solve_flag] = obj.solve_with_yalmip(constraints_all, L, obj.constParams.yalmip_options);

            % Get the values of the variables
            obj.this_iter.vars        = obj.assign_values(obj.scp_prob.sdp_vars);
            obj.this_iter.J0          = value(J0);
            obj.this_iter.P_cvx       = value(P_cvx);
            obj.this_iter.L           = value(L);
            obj.this_iter.DeltaZ      = value(DeltaZ);
            obj.this_iter.trust_region_constraint_lhs = value(trust_region_constraint_lhs);
            obj.this_iter.slack_noncvx_eq      = value(obj.scp_prob.slack_noncvx_eq);
            obj.this_iter.slack_noncvx_ineq    = value(obj.scp_prob.slack_noncvx_ineq);
            obj.this_iter.time_solver = diagnostics.solvertime;
            obj.this_iter.time_yalmip = diagnostics.yalmiptime;
        end

        % Compute the augmented penalty function
        function P = penalty_function(obj, xi, zeta, options)

            arguments
                obj
                xi
                zeta
                options.DeltaZ = []
                options.trust_region_constraint_lhs = []
            end

            % Extract for readability
            lambda  = obj.this_iter.lambda;
            mu      = obj.this_iter.mu;
            w       = obj.this_iter.w;
            w_const = obj.constParams.w_const;
            % w_TR    = obj.constParams.w_TR;
            
            zeta_pos = obj.pos(zeta);

            switch obj.constParams.penalty_method
                case {'AL'}
                    P = dot(lambda, xi) + w/2 * dot(xi,xi) + dot(mu,zeta_pos) + w/2 * dot(zeta_pos, zeta_pos);

                % l1 penalty as seen in Mao et al. 2018
                case {'L1'}
                    P = w_const * (norm(xi,1) + sum(zeta_pos));

                % Add AL and l1 penalty; take sqrt(w) to ensure that the penalty is scaled correctly
                case {'ALwithL1'}
                    P = dot(lambda, xi) + w/2 * dot(xi,xi) + dot(mu, zeta_pos) + w/2 * dot(zeta_pos, zeta_pos) + ...
                        sqrt(w) * (norm(xi,1) + sum(zeta_pos));

                % Additional soft trust region penalty
                case {'AL_with_softTR'}
                    P = dot(lambda, xi) + w/2 *dot(xi,xi) + dot(mu, zeta_pos) + w/2 * dot(zeta_pos, zeta_pos) ...
                        + sqrt(w) * obj.pos(options.trust_region_constraint_lhs);
            end

        end

        % Postprocessing of results after each subproblem solve
        function [next_iter] = postprocess(obj, verbose, clear_every_iter)

            % In case of numerical issues, reduce the penalty weight and resolve
            if obj.this_iter.solve_flag == 4 && strcmp(obj.constParams.numerical_issue_handling, 'resolve_with_reduced_penalty')
                next_iter = SCPIter(obj.this_iter.iter + 1, obj);
                next_iter.w = next_iter.w / 2;
                disp("Going to next iteration with reduced penalty w = " + next_iter.w);
                return
            end
            
            % Get the nonconvex constraint violations
            obj.this_iter.g = obj.scp_prob.noncvx_eq(obj.this_iter.vars);
            obj.this_iter.h = obj.scp_prob.noncvx_ineq(obj.this_iter.vars);

            % Stack the constraint violations 
            obj.this_iter.chi = norm([obj.this_iter.g; max(0, obj.this_iter.h)]);

            % Augmented cost of the newest solution evaluated with the nonconvex constraints
            obj.this_iter.P = obj.penalty_function(obj.this_iter.g, obj.this_iter.h, DeltaZ=obj.this_iter.DeltaZ);
            J = obj.this_iter.J0 + obj.this_iter.P;

            % Augmented cost of the reference with the current w, lambda, mu values
            obj.this_iter.P_ref = obj.penalty_function(obj.this_iter.g_ref, obj.this_iter.h_ref, DeltaZ=obj.this_iter.DeltaZ);
            J_ref = obj.this_iter.J0_ref + obj.this_iter.P_ref;

            % Evaluate the nonconvex/relaxed versions of improvement in augmented cost
            obj.this_iter.delJ = J_ref - J;
            obj.this_iter.delL = J_ref - obj.this_iter.L;

            % Initialize the next iteration
            next_iter = SCPIter(obj.this_iter.iter + 1, obj);

            % DeltaL should be positive by asssumption when the constraints are linearized exactly
            if obj.this_iter.delL < 0 && strcmp(obj.constParams.linearization, 'exact')
                obj.disp_verbose(verbose, ...
                "DeltaL = " + obj.this_iter.delL + " < 0. Constraint linearization might be incorrect. Proceeding...");
                % keyboard
            end

            if obj.this_iter.delL == 0
                obj.this_iter.rho = 1;
                obj.fprintf_verbose(verbose, '\nDeltaL = 0\n')
                % keyboard
            else 
                obj.this_iter.rho = obj.this_iter.delJ / obj.this_iter.delL; % ratio of actual cost reduction to predicted reduction
            end

            next_iter = obj.update_trust_region(next_iter);

            if (strcmp(obj.constParams.linearization, 'exact') && obj.this_iter.rho >= obj.constParams.rho0) || ...
                (strcmp(obj.constParams.linearization, 'inexact') && 1 - obj.constParams.eta0 <= obj.this_iter.rho && obj.this_iter.rho <= 1 + obj.constParams.eta0)
                % step is accepted
                obj.this_iter.update_ref = true;
                next_iter.ref_vars = obj.this_iter.vars;
                next_iter.g_ref    = obj.this_iter.g;
                next_iter.h_ref    = obj.this_iter.h;
                next_iter.J0_ref   = obj.this_iter.J0;

                if ~obj.use_optimizer
                    obj.scp_prob.update_parameters(next_iter.ref_vars);
                end
                
                % (Eq 15)
                if abs(obj.this_iter.delJ) < min(obj.this_iter.delta, obj.this_iter.eta * obj.this_iter.chi) 

                    obj.this_iter.update_multipliers = true;

                    if obj.constParams.superlinear && obj.this_iter.eta == Inf
                        next_iter.eta = abs(obj.this_iter.delJ) / obj.this_iter.chi;
                    end
                    
                    next_iter = obj.update_multipliers(obj.this_iter.g, obj.this_iter.h, obj.this_iter.delJ, next_iter);
                end

            end

            if ~obj.use_optimizer && clear_every_iter
                yalmip('clear');
                obj.scp_prob.set_fixed_sdp_objects();
                obj.scp_prob.update_constraints(next_iter.ref_vars);
                return
            end

            if obj.this_iter.update_ref
                obj.scp_prob.update_constraints(next_iter.ref_vars);
            end
        end

        % Update multipliers as well as weights after each successful iteration
        function next_iter = update_multipliers(obj, g, h, delJ, next_iter)

            % multiplier update Eqn (3)
            next_iter.lambda = obj.this_iter.lambda + obj.this_iter.w * g;
            next_iter.mu     = max(0, obj.this_iter.mu + obj.this_iter.w * h);
            next_iter.w      = min(obj.constParams.beta * obj.this_iter.w, obj.constParams.w_max);

            % stationarity tolerance update Eqn (16)
            if obj.this_iter.delta == Inf
                next_iter.delta = abs(delJ);
            else
                next_iter.delta = obj.constParams.gamma * obj.this_iter.delta;
            end
        end

        % Update the trust region radius based on rho; Eqn (17)
        function next_iter = update_trust_region(obj, next_iter)

            switch obj.constParams.linearization
                case 'exact'
                    if obj.this_iter.rho < obj.constParams.rho1 % shrink trust region
                        next_iter.r = max(obj.this_iter.r / obj.constParams.alpha1, obj.constParams.r_min);
                    elseif obj.this_iter.rho >= obj.constParams.rho2 % expand trust region
                        next_iter.r = min(obj.this_iter.r * obj.constParams.alpha2, obj.constParams.r_max);
                    end
                case 'inexact'
                    if 1 - obj.constParams.eta2 <= obj.this_iter.rho && obj.this_iter.rho <= 1 + obj.constParams.eta2
                        next_iter.r = min(obj.this_iter.r * obj.constParams.alpha2, obj.constParams.r_max); % expand trust region
                    elseif 1 - obj.constParams.eta1 <= obj.this_iter.rho && obj.this_iter.rho <= 1 + obj.constParams.eta1
                        % do nothing
                    else
                        next_iter.r = max(obj.this_iter.r / obj.constParams.alpha1, obj.constParams.r_min); % shrink trust region
                    end
            end

        end

        function is_stalled = check_stall(obj, iter)
            if iter > 5
                delJ_hist = [obj.iter_hist.delJ];
                chi_hist = [obj.iter_hist.chi];
                delJ_last5 = delJ_hist(end-4:end);
                chi_last5 = chi_hist(end-4:end);
                % If the last 5 iterations have similar delJ and chi, then the algorithm is stalled
                is_stalled = all(abs(diff(delJ_last5)) < 1E-6) && all(abs(diff(chi_last5)) < 1E-6);
            else
                is_stalled = false;
            end
        end

        function has_numerical_issues = check_numerical(obj, iter)
            if iter > 5
                % if more than three out of the last five iterations have solve_flag = 4, then there are numerical issues
                solve_flag_hist = [obj.iter_hist.solve_flag];
                has_numerical_issues = sum(solve_flag_hist(end-4:end) == 4) > 3;
            else
                has_numerical_issues = false;
            end
        end

        function is_close_to_optimum = check_close_to_optimum(obj, iter, last_checked_close_to_optimum)
            if last_checked_close_to_optimum + 10 < iter
                % if more than eight out of the last 10 iterations satisfy the feasibility tolerance
                % and the objective function isn't decreasing very smoothly, then the algorithm is close to the optimum
                chi_hist = [obj.iter_hist.chi];
                delJ_hist = [obj.iter_hist.delJ];
                is_close_to_optimum = sum(chi_hist(end-9:end) < obj.constParams.tol_feas) > 8;
                is_close_to_optimum = is_close_to_optimum && sum(delJ_hist(end-9:end) > 0) > 9;
            else
                is_close_to_optimum = false;
            end
        end

        % Check stopping criteria for the iteration; Eqn (18)
        function stop = check_stopping_criteria(obj)
            stop = (abs(obj.this_iter.delJ) <= obj.constParams.tol_opt) ...
                    && (obj.this_iter.chi <= obj.constParams.tol_feas) ...
                    && (norm(obj.this_iter.DeltaZ, Inf) <= obj.constParams.tol_change);
        end
                
        % Create a report at the end of the optimization process
        function create_report(obj, iter, converged)
            obj.report.iters = iter;
            obj.report.time = datetime('now') - obj.starttime;
            obj.report.converged = converged;
        end

        function constraintLHS = get_trust_region_constraint_lhs(obj, vars)
            constraintLHS = [];
            for i = 1:length(obj.var_names)
                field = obj.var_names{i};
                if obj.scp_prob.impose_trust_region_struct.(field)
                        constraintLHS = [constraintLHS
                            obj.scp_prob.trust_region_scaling.(field) ...
                            .* norm(reshape(vars.(field) - obj.this_iter.ref_vars.(field), [], 1), obj.constParams.trust_region_norm) ...
                           - obj.this_iter.r
                        ];
                    
                end
            end
        end

        % Get the difference in the variables from the reference as a concatenated vector
        function DeltaZ = getDeltaZ(obj, vars)

            DeltaZ = [];
            for i = 1:length(obj.var_names)
                field = obj.var_names{i};
                if obj.scp_prob.impose_trust_region_struct.(field)
                    DeltaZ = [DeltaZ; reshape(vars.(field) - obj.this_iter.ref_vars.(field), [], 1)];
                end
            end
        end
           
        function plot_iter_history(obj)
            figure;
            tiledlayout(2,1);
            
            nexttile;
            hold on;
            abs_delJ = abs([obj.iter_hist.delJ]);
            plot(abs_delJ, 'k', 'HandleVisibility', 'off');

            % Reference updates
            delJ_ref_updates = abs_delJ;
            delJ_ref_updates(~[obj.iter_hist.update_ref]) = NaN;
            semilogy(delJ_ref_updates, 'rx', 'LineStyle', 'none', 'DisplayName', 'Reference Update');
            % Multiplier updates
            delJ_multiplier_updates = abs_delJ;
            delJ_multiplier_updates(~[obj.iter_hist.update_multipliers]) = NaN;
            semilogy(delJ_multiplier_updates, 'bo', 'LineStyle', 'none', 'DisplayName', 'Multiplier Update');
            semilogy([obj.iter_hist.delta], 'DisplayName', '$\delta$')
            if obj.constParams.superlinear
                semilogy([obj.iter_hist.eta] .* [obj.iter_hist.chi], 'DisplayName','$\eta \cdot\chi$')
            end

            yline(obj.constParams.tol_opt, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
            xlabel('Iteration')
            ylabel('$|\Delta J|$')
            legend()
            yscale log

            nexttile;
            hold on;
            semilogy([obj.iter_hist.chi], 'k', 'HandleVisibility', 'off');

            % Reference updates
            chi_ref_updates = [obj.iter_hist.chi];
            chi_ref_updates(~[obj.iter_hist.update_ref]) = NaN;
            semilogy(chi_ref_updates, 'rx', 'LineStyle', 'none', 'DisplayName', 'Reference Update');

            % Multiplier updates
            chi_multiplier_updates = [obj.iter_hist.chi];
            chi_multiplier_updates(~[obj.iter_hist.update_multipliers]) = NaN;
            semilogy(chi_multiplier_updates, 'bo', 'LineStyle', 'none', 'DisplayName', 'Multiplier Update');
            
            yline(obj.constParams.tol_feas, 'k--', 'LineWidth', 1.5, 'HandleVisibility', 'off');
            xlabel('Iteration')
            ylabel('$\chi$')

            yscale log
        end

        function print_iteration_header(obj, verbose)
            obj.fprintf_verbose(verbose, '\nIteration       delJ         chi          J0         P         Ref. Upd.      Mult. Upd.\n');
        end

        function print_iteration_info(obj, verbose)
            header_frequency = 10;

            if mod(obj.this_iter.iter-1, header_frequency) == 0
                obj.print_iteration_header(verbose);
            end
            if obj.this_iter.update_ref
                update_ref_yn = 'Y';
            else
                update_ref_yn = 'N';
            end
            if obj.this_iter.update_multipliers
                update_multipliers_yn = 'Y';
            else
                update_multipliers_yn = 'N';
            end
 
            obj.fprintf_verbose(verbose, '%8d ', obj.this_iter.iter);
            if abs(obj.this_iter.delJ) < obj.constParams.tol_opt
                obj.cprintf('key', '%12.4e ', obj.this_iter.delJ);
            else
                obj.cprintf('text', '%12.4e ', obj.this_iter.delJ);
            end
            if obj.this_iter.chi < obj.constParams.tol_feas
               obj.cprintf('key', '%12.4e ', obj.this_iter.chi);
            else
                obj.cprintf('text', '%12.4e ', obj.this_iter.chi);
            end
            obj.cprintf('text', '%12.4e ', obj.this_iter.J0);
            obj.cprintf('text', '%12.4e ', obj.this_iter.P);
            obj.cprintf('text', '      %s            %s\n', update_ref_yn, update_multipliers_yn);

        end

        function plot_time(obj)
            % For each iteration, plot a bar graph of the time spent in each phase
            time_yalmip = [obj.iter_hist.time_yalmip];
            time_solver = [obj.iter_hist.time_solver];
            time_postprocessing = [obj.iter_hist.time_postprocessing];

            figure;
            bar([time_yalmip; time_solver; time_postprocessing]', 'stacked');
            xlabel('Iteration')
            ylabel('Time (s)')
            legend('YALMIP', 'Solver', 'Postprocessing')
        end

        function print_time_info(obj)
            time_yalmip = [obj.iter_hist.time_yalmip];
            time_solver = [obj.iter_hist.time_solver];
            time_postprocessing = [obj.iter_hist.time_postprocessing];

            time_yalmip_avg = mean(time_yalmip);
            time_solver_avg = mean(time_solver);
            time_postprocessing_avg = mean(time_postprocessing);

            fprintf('Average time spent in YALMIP: %f seconds\n', time_yalmip_avg);
            fprintf('Average time spent in solver: %f seconds\n', time_solver_avg);
            fprintf('Average time spent in postprocessing: %f seconds\n', time_postprocessing_avg);
        end

        function plot_trust_region_history(obj)
            figure
            plot([obj.iter_hist.r]);
            xlabel('Iteration')
            ylabel('Trust Region Radius $r$')
            yscale log
        end

        function plot_objective_history(obj)
            figure
            plot([obj.iter_hist.J0]);
            xlabel('Iteration')
            ylabel('Objective Function $J_0$')
        end

        function plot_augmented_cost_history(obj)
            figure
            plot([obj.iter_hist.L]);
            xlabel('Iteration')
            ylabel('Augmented Cost $L$')
        end

        function plot_penalty_history(obj)
            figure
            plot([obj.iter_hist.P]);
            xlabel('Iteration')
            ylabel('Penalty $P$')
        end

        function vars_double = assign_values(obj, vars)
            % Assign values to sdpvar objects after solving the optimization problem
            for i = 1:length(obj.var_names)
                field = obj.var_names{i};
                vars_double.(field) = value(vars.(field));
            end
        end

    end

    methods (Static)

        % Helper function to ensure that the output is positive
        function out = pos(x)
            % sdpvar is constrained to be positive, so we can just return x
            if isa(x, 'sdpvar')
                out = x;
            elseif isa(x, 'double')
                out = max(0, x);
            end
        end

        function fprintf_verbose(verbose, varargin)
            if verbose
                fprintf(varargin{:});
            end
        end

        function disp_verbose(verbose, varargin)
            if verbose
                disp(varargin{:});
            end
        end

        function [diagnostics, solve_flag] = solve_with_yalmip(constraints, objective, yalmip_options)
            arguments
                constraints
                objective
                yalmip_options = sdpsettings();
            end
            diagnostics = optimize(constraints, objective, yalmip_options);
            solve_flag = diagnostics.problem;
            switch solve_flag
                case 0 % Solved

                case 1 % Infeasible. resolve problem without objective 
                    fprintf('not solved: %s\n', diagnostics.info)
                    feasibility_problem = optimize(constraints, [], yalmip_options);
                    if feasibility_problem.problem == 0
                        fprintf('Feasibility problem solved. There may be a problem in the objective. \n')                        
                    else
                        fprintf('Feasibility problem not solved. The problem is infeasible. \n')
                    end
                case 4 % Numerical problems
                    fprintf('Numerical problems: MOSEK Error Code %i\n', diagnostics.solveroutput.r)
                otherwise
                    fprintf('not solved: %s\n', diagnostics.info)
            end
            % if yalmip_options.savesolveroutput && ~isempty(diagnostics.solveroutput.res.rmsg)
                % fprintf(diagnostics.solveroutput.res.rmsg)
            % end

        end

        count = cprintf(style,format,varargin) % Print in color. Defined in separate file
        
    end
        
end
