%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ExampleClass1Optimizer.m is a class that inherits from ExampleClass1.m and
% uses the optimizer functionality of yalmip for acclerated computation.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef ExampleClass1Optimizer < ExampleClass1

    properties
        changing_params = struct();
        is_used_in_constraint_update = struct('z', true);
    end
    
    methods
        % Constructor
        function obj = ExampleClass1Optimizer()
            obj@ExampleClass1();
        end

        function val = get_use_optimizer(obj)
            val = true;
        end

        function set_changing_parameters_sdp(obj)
            obj.changing_params.f_ref = sdpvar(1);
            obj.changing_params.dfdz_ref = sdpvar(2,1);
        end

        function p = get_changing_parameters(obj, ref_vars)
            p.f_ref = obj.f(ref_vars.z);
            p.dfdz_ref = obj.dfdz(ref_vars.z);
        end
            
        function constraints = get_changing_constraints(obj, vars, ref_vars)
            obj.set_changing_parameters_sdp()
            p = obj.changing_params;
            constraints = [
                p.f_ref + p.dfdz_ref'*(vars.z - ref_vars.z) == obj.slack_noncvx_eq
            ];
        end

    end
end
