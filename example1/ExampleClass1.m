%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ExampleClass1.m Implemention of Example 1 from Oguri 2023.
% This file should be used as a template for creating new SCP problems.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef ExampleClass1 < SCPProblem

    properties
        % impose_trust_region_struct, init_guess_struct, and SCvxParams must be specified
        % For each variable, specify whether to impose a trust region
        impose_trust_region_struct = struct('z', true);
        % For each variable, specify an initial guess value
        init_guess_struct = struct('z', [1.5; 1.5]);
        
        % Utility properties (not mandatory)
        f = @(z) z(2) - z(1)^4 - 2*z(1)^3 + 1.2*z(1)^2 + 2*z(1);
        dfdz = @(z) [-4*z(1)^3 - 6*z(1)^2 + 2.4*z(1) + 2; 1];
        fig
    end

    properties
        f_ref
        dfdz_ref
    end
    
    methods
        % Constructor
        function obj = ExampleClass1()
            % Must inherit from SCPProblem class
            obj@SCPProblem();
            % Must initialize the parent class
            obj.initialize();
        end

        % Define optimization variables
        function vars = define_vars(obj)
            vars.z = sdpvar(2,1); % Define z as a 2x1 optimization variable
        end

        % Objective function
        function J0 = objective(obj, vars)
            J0 = sum(vars.z);
        end

        % Convex equality constraints
        function constraints = convex_eq(obj, vars)
            constraints = [];
        end

        % Convex inequality constraints
        function constraints = convex_ineq(obj, vars)
            constraints = [
                -2 <= vars.z
                vars.z <= 2
                -vars.z(2) - 4/3 * vars.z(1) - 2/3 <= 0
            ];
        end

        function update_parameters(obj, ref_vars)
            obj.f_ref = obj.f(ref_vars.z);
            obj.dfdz_ref = obj.dfdz(ref_vars.z);
        end
                
        % Non-convex equality constraints
        function constraintLHS = noncvx_eq(obj, vars)
            constraintLHS = obj.f(vars.z);
        end

        % Non-convex equality constraints (relaxed)
        function constraintLHS = noncvx_eq_relaxed(obj, vars, ref_vars)
            constraintLHS = obj.f_ref + obj.dfdz_ref.'*(vars.z - ref_vars.z);
        end
        
        % Non-convex inequality constraints
        function constraintLHS = noncvx_ineq(obj, vars)
            constraintLHS = [];
        end

        % Non-convex inequality constraints (relaxed)
        function constraintLHS = noncvx_ineq_relaxed(obj, vars, ref_vars)
            constraintLHS = [];
        end

        % Post-iteration code here (optional)
        function post_iteration(obj)
            if obj.scp.this_iter.iter == 1
                obj.fig = figure;
                obj.plot_setting()
            end
            obj.fig;
            z = obj.scp.this_iter.vars.z;
            plot(z(1), z(2), 'b*')
            hold on
            drawnow
            % exportgraphics(gcf,'example1.gif','Append',true);
        end

        function plot_setting(obj)
            z1_range = -2:0.1:2;
            cnst_b = @(z1) z1.^4 + 2*z1.^3-1.2*z1.^2-2*z1;
            cnst_c = @(z1) -4/3*z1-2/3;   
            hold on
            plot(z1_range, cnst_b(z1_range), 'k')
            plot(z1_range, cnst_c(z1_range), 'k')
            x = [z1_range, fliplr(z1_range)];
            y = [-2*ones(1,length(z1_range)), fliplr(cnst_c(z1_range))];
            fill(x, y, 'r', 'FaceAlpha', 0.3)
    
            xlabel('$z_1$')
            ylabel('$z_2$')
            xlim([-2, 2])
            ylim([-2, 2])
        end

        function plot_contour(obj)
            z1_range = -2:0.1:2;
            z2_range = -2:0.1:2;
            [Z1, Z2] = meshgrid(z1_range, z2_range);
            J0 = Z1 + Z2;
            contourf(Z1, Z2, J0, 100, 'LineStyle', 'none')
            cb = colorbar;
            cb.Label.String = '$J(z)$';
            cb.Label.Interpreter = 'latex';
            hold on
            obj.plot_setting()
        end
    end
end
