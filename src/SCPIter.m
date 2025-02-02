classdef SCPIter < handle
    properties
        iter {mustBeInteger}

        % Set in the previous iteration
        lambda {isreal}
        mu {isreal}
        w {isreal}
        delta
        r
        eta

        ref_vars struct
        J0_ref
        g_ref
        h_ref

        % Set after solving the subproblem
        vars        = struct()
        J0          = NaN
        P_cvx       = NaN
        L           = NaN
        DeltaZ      = NaN
        trust_region_constraint_lhs = NaN
        slack_noncvx_eq = NaN
        slack_noncvx_ineq = NaN
        solve_flag  = NaN
        time_yalmip = NaN
        time_solver = NaN

        % Set in the post-processing step
        P_ref               = NaN
        g                   = NaN
        h                   = NaN
        P                   = NaN
        delJ                = NaN
        delL                = NaN
        rho                 = NaN
        chi                 = NaN
        update_ref          = false
        update_multipliers  = false
        time_postprocessing = NaN

    end

    methods
        function obj = SCPIter(iter, scp)
            arguments
                iter {mustBeInteger}
                scp SCvxStar
            end
            obj.iter = iter;
            if obj.iter > 1

                obj.lambda   = scp.this_iter.lambda;
                obj.mu       = scp.this_iter.mu;
                obj.w        = scp.this_iter.w;
                obj.delta    = scp.this_iter.delta;
                obj.r        = scp.this_iter.r;
                obj.eta      = scp.this_iter.eta;

                obj.ref_vars = scp.this_iter.ref_vars;
                obj.J0_ref   = scp.this_iter.J0_ref;
                obj.g_ref    = scp.this_iter.g_ref;
                obj.h_ref    = scp.this_iter.h_ref;

                return
            end

            obj.w = scp.constParams.w_init;
            obj.r = scp.constParams.r_init;
            obj.delta = Inf;
            if scp.scp_prob.size_noncvx_eq > 0
                obj.lambda = zeros(scp.scp_prob.size_noncvx_eq, 1);
            else
                obj.lambda = [];
            end
            if scp.scp_prob.size_noncvx_ineq > 0
                obj.mu = zeros(scp.scp_prob.size_noncvx_ineq, 1);
            else
                obj.mu = [];
            end
            obj.eta = Inf;

        end
        
    end
end