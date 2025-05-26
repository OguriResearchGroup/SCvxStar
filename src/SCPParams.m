classdef SCPParams < handle

    properties
        r_init {mustBePositive}     = 0.1; % initial trust region radius
        w_init {mustBePositive}     = 100; % initial penalty parameter
        tol_opt {mustBePositive}    = 1e-5; % optimality tolerance
        tol_feas {mustBePositive}   = 1e-5; % feasibility tolerance
        tol_change {mustBePositive} = Inf; % variable change tolerance
        rho0                        = 0; % minimum value of rho to accept solution and update the reference (for exact linearization)
        rho1                        = 0.25; % rho0 < rho < rho1: decrease trust region radius (for exact linearization)
        rho2                        = 0.7; % rho2 < rho: increase trust region radius (for exact linearization)
        eta0                        = 1; % accept solution and update the reference if rho is within 1 - eta0 and 1 + eta0 (for inexact linearization)
        eta1                        = 0.75; % decrease trust region radius if rho is not within 1 - eta1 and 1 + eta1 (for inexact linearization)
        eta2                        = 0.1; % eta2 < eta1: increase trust region radius if rho is within 1 - eta2 and 1 + eta2 (for inexact linearization)
        alpha1 {mustBePositive}     = 2; % divide trust region radius by alpha1
        alpha2 {mustBePositive}     = 3; % multiply trust region radius by alpha2
        beta {mustBePositive}       = 2; % multiply penalty parameter by beta
        gamma {mustBePositive}      = 0.9; % stationarity parameter
        r_min {mustBePositive}      = 1e-6; % minimum trust region radius
        r_max {mustBePositive}      = 1; % maximum trust region radius
        w_max {mustBePositive}      = 1e8; % maximum penalty parameter
        k_max {mustBeInteger}       = 2000; % maximum number of iterations
        w_const {mustBePositive}    = 1e4; % penalty when using constant penalty as in Mao et al.
        w_TR {mustBePositive}       = 100; % penalty for soft trust region. not used if trust region is not imposed on any variable
        penalty_method {ismember(penalty_method, {'AL', 'ALwithL1', 'L1', 'AL_with_softTR'})} = 'AL';
        superlinear {mustBeNumericOrLogical} = false; % superlinear convergence
        yalmip_options struct                = sdpsettings('verbose', 0, 'solver', 'mosek', 'warning', 1, 'beeponproblem', 1, 'savesolveroutput', 1);
        soft_TR {mustBeNumericOrLogical}     = false; % use soft trust region, activated only when the original two convergence criteria are met
        trust_region_norm                    = inf; % norm used for trust region
        linearization {ismember(linearization, {'exact', 'inexact'})} = 'exact'; % linearization method
        % how to handle numerical issues; accept or resolve with smaller penalty
        numerical_issue_handling {ismember(numerical_issue_handling, {'accept', 'resolve_with_reduced_penalty'})} = 'accept'; 
    end

    methods
        function obj = SCPParams()
        end

        function set.rho0(obj, value)
            assert(value >= 0, 'rho0 must be non-negative');
            assert(value < obj.rho1, 'rho0 must be less than rho1');
            obj.rho0 = value;
        end

        function set.rho1(obj, value)
            assert(value > obj.rho0, 'rho1 must be greater than rho0');
            assert(value < obj.rho2, 'rho1 must be less than rho2');
            obj.rho1 = value;
        end

        function set.rho2(obj, value)
            % assert(value > obj.rho1, 'rho2 must be greater than rho1');
            obj.rho2 = value;
        end

        function set.eta0(obj, value)
            assert(1 >= value && value > 0, 'eta0 must be 1 >= eta0 > 0');
            % assert(value > obj.eta1, 'eta0 must be greater than eta1');
            obj.eta0 = value;
        end

        function set.eta1(obj, value)
            assert(1 > value && value > 0, 'eta1 must be between 0 and 1');
            % assert(obj.eta0 > value && value > eta2, 'eta1 must be between eta0 and eta2');
            obj.eta1 = value;
        end

        function set.eta2(obj, value)
            assert(1 > value && value > 0, 'eta2 must be between 0 and 1');
            % assert(eta1 > value, 'eta2 must be less than eta1');
            obj.eta2 = value;
        end
    end
end