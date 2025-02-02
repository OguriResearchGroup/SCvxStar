%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SCPPProblem
% This defines a template class with basic methods for handling constraints
% The user is not supposed to call this class directly, but rather
% through defining their own class that inherits this class
% See ExampleClass1.m for an example of how to define a new SCP problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef SCPProblem < handle

    % These properties will be set by the user
    properties (Abstract)
        impose_trust_region_struct {is_struct_with_bool}
        init_guess_struct {is_struct_with_double}
    end

    % These properties are internally set by the class
    properties
        sdp_vars {is_struct_with_sdpvar}
        sdp_objective {is_sdpvar}
        sdp_convex_eq {is_sdpvar}
        sdp_convex_ineq {is_sdpvar}
        sdp_noncvx_eq_relaxed {is_sdpvar}
        sdp_noncvx_ineq_relaxed {is_sdpvar}
        sdp_convexified_exact {is_sdpvar}
        size_noncvx_eq double
        size_noncvx_ineq double
        slack_noncvx_eq = []; % yalmip objects or empty
        slack_noncvx_ineq = []; % yalmip objects or empty
        trust_region_scaling = []; % Trust region scaling for each field. Can be scalar or vector of the same size as the number of variables
            % trust_region_scaling can be updated iteratively if the trust region is adaptive (e.g. state-dependent)
        scp SCvxStar % SCvxStar object
        sol % SCvxStar solution
    end

    properties (Dependent)
        use_optimizer
    end

    methods

        % Constructor
        function obj = SCPProblem()

        end

        function initialize(obj)
            % Define optimization variables, objective, and constraints as yalmip objects
            yalmip('clear');

            % Automatically determine the size of nonconvex constraints
            obj.size_noncvx_eq = length(obj.noncvx_eq(obj.init_guess_struct));
            obj.size_noncvx_ineq = length(obj.noncvx_ineq(obj.init_guess_struct));

            obj.set_fixed_sdp_objects();

            if ~obj.use_optimizer
                obj.update_parameters(obj.init_guess_struct);
                obj.update_constraints(obj.init_guess_struct);
            end

            if isempty(obj.trust_region_scaling) % If trust region scaling is not provided, set all fields to 1
                obj.trust_region_scaling = struct();
                fields = fieldnames(obj.sdp_vars);
                for i = 1:length(fields)
                    % Get the field name
                    field = fields{i};

                    if obj.impose_trust_region_struct.(field)
                        obj.trust_region_scaling.(field) = 1;
                    end
                end
            else % If trust region scaling is provided, check that all fields are present and the values are either
                % scalar or vector of the same size as the number of variables
                fields = fieldnames(obj.sdp_vars);
                for i = 1:length(fields)
                    % Get the field name
                    field = fields{i};

                    if ~obj.impose_trust_region_struct.(field)
                        continue;
                    end

                    if ~isfield(obj.trust_region_scaling, field)
                        error('Trust region scaling not provided for field %s', field);
                    end

                    if ~isnumeric(obj.trust_region_scaling.(field)) || ...
                        (~isscalar(obj.trust_region_scaling.(field)) && all(size(obj.trust_region_scaling.(field)) ~= size(obj.sdp_vars.(field))))
                        error('Trust region scaling for field %s must be a scalar or vector of the same size as the number of variables', field);
                    end
                end
                
            end
        end

        function set_fixed_sdp_objects(obj)
            % Define slack variables for nonconvex constraints
            if obj.size_noncvx_eq > 0
                obj.slack_noncvx_eq = sdpvar(obj.size_noncvx_eq, 1);
            end

            if obj.size_noncvx_ineq > 0
                obj.slack_noncvx_ineq = sdpvar(obj.size_noncvx_ineq, 1);
            end

            obj.sdp_vars        = obj.define_vars();
            obj.sdp_objective   = obj.objective(obj.sdp_vars);
            obj.sdp_convex_eq   = obj.convex_eq(obj.sdp_vars);
            obj.sdp_convex_ineq = obj.convex_ineq(obj.sdp_vars);
        end

        function varargout = solve(obj, options)
            arguments
                obj
                options.scp_params = SCPParams()
                options.verbose = true
                options.save_bool = true; % If converged, save results to .mat file in the logs folder
                options.start_from_constructor = true; % If true, initialize the problem from the constructor
                options.clear_every_iter = false; % If true, clear the YALMIP variables every iteration
            end

            if options.start_from_constructor
                % Initialize the SCvxStar object
                obj.scp = SCvxStar(obj, options.scp_params, obj.use_optimizer);
            end

            % Solve the problem
            [converged] = obj.scp.solve(verbose = options.verbose, clear_every_iter = options.clear_every_iter);

            obj.sol = obj.scp.this_iter.vars;
            
            if converged && options.save_bool
                obj.save_to_mat();
            end

            if nargout > 0
                varargout{1} = converged;
            end
        end

        function update_constraints(obj, ref_vars, options)
            arguments
                obj
                ref_vars
                options.update_parameters = false
            end
            % Set the nonconvex constraints using the reference variables
            % This is called when the reference variables change
            % Inherit this function and modify to your needs; e.g. calculating STMs
            obj.sdp_noncvx_eq_relaxed = obj.noncvx_eq_relaxed(obj.sdp_vars, ref_vars) == obj.slack_noncvx_eq;
            obj.sdp_noncvx_ineq_relaxed = obj.noncvx_ineq_relaxed(obj.sdp_vars, ref_vars) <= obj.slack_noncvx_ineq;
            obj.sdp_convexified_exact = obj.convexified_exact(obj.sdp_vars, ref_vars);
        end

        function save_to_mat(obj)
            % Save the object as a .mat file
            obj.sol.objective = value(obj.sdp_objective);
            yalmip('clear'); % YALMIP variables cannot be saved
            filename = string(datetime("now", 'Format', 'yyyy-MM-dd_HH:mm')) + ".mat";
            full_path = fullfile(fileparts(mfilename('fullpath')), "..", "logs", filename);

            % If logs directory does not exist, create it
            if ~exist(fullfile(fileparts(mfilename('fullpath')), "..", "logs"), 'dir')
                mkdir(fullfile(fileparts(mfilename('fullpath')), "..", "logs"));
            end

            save(full_path, 'obj');
            disp('Saved to SCvxStar/logs/' + filename);
            disp('To load the object, use the "load_scp" function');
        end

        function val = get.use_optimizer(obj)
            val = get_use_optimizer(obj);
        end

        function val = get_use_optimizer(obj)
            val = false;
        end
    end

    % These properties are set by the user
    methods (Abstract)

        % Define optimization variables
        vars = define_vars(obj);

        % Objective function
        J0 = objective(obj, vars);

        % Convex equality constraints
        constraints = convex_eq(obj, vars);

        % Convex inequality constraints
        constraints = convex_ineq(obj, vars);

        % Non-convex equality constraints
        constraintLHS = noncvx_eq(obj, vars);

        % Non-convex inequality constraints
        constraintLHS = noncvx_ineq(obj, vars);

        % Linearized non-convex equality constraints
        constraintLHS = noncvx_eq_relaxed(obj, vars, ref_vars);

        % Linearized non-convex inequality constraints
        constraintLHS = noncvx_ineq_relaxed(obj, vars, ref_vars)

    end

    % Optional methods that can be overwritten by the user
    methods
        function constraint = convexified_exact(obj, vars, ref_vars)
            % Convex/convexified constraints that are imposed exactly but change with the reference variables
            % WARNING: In general, this breaks the convergence guarantee of SCvx/SCvx* and causes Delta L to become negative
            % Use with caution
            constraint = [];
        end

        function pre_iteration(obj)
            % Do nothing by default
        end

        % Things to do after each the postprocess step of each iteration (e.g. plotting)
        function post_iteration(obj)
            % Do nothing by default
        end

        %% The following functions are used when using the optimizer (i.e. use_optimizer = true)
        function set_changing_parameters_sdp(obj)
            % Set the changing parameters as sdpvars
        end

        function p = get_changing_parameters(obj, ref_vars)
            % Get the changing parameters
            p = [];
        end

        function constraints = get_changing_constraints(obj, vars, ref_vars)
            % Get the constraints that depend on the changing parameters
            constraints = [];
        end

        function [p, constraints] = get_parameterized_constraints(obj, vars)
            % Get the problem parameters and the constraints that depend on the problem parameters
            p = [];
            constraints = [];
        end

        function update_parameters(obj, ref_vars)
            % Update the parameters that depend on the reference variables
        end

    end

end

%% Helper functions for type checking
function out = is_sdpvar(x)
    out = isa(x, 'sdpvar');
end

function out = is_struct_with_sdpvar(x)
    out = isstruct(x) && all(structfun(@is_sdpvar, x));
end

function out = is_struct_with_bool(x)
    out = isstruct(x) && all(structfun(@islogical, x));
end

function out = is_struct_with_double(x)
    out = isstruct(x) && all(structfun(@isnumeric, x));
end