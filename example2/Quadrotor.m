%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Quadrotor class definition
% Requires the DynamicalSystem class 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef Quadrotor < astro.DynamicalSystem

    properties (SetAccess=immutable)
        m % mass
        k_D % drag coefficient
        g % gravity vector
    end
    properties (Constant)
        is_control_affine = true;
        nx = 6;
        nu = 3;
        B_impulse = [];
    end

    methods
        function obj = Quadrotor(m, k_D, g)
            obj.m = m;
            obj.k_D = k_D;
            obj.g = g;
        end   

        function dxdt = EOM(obj, t, x)
            v = x(4:6);
            dxdt = [v; - obj.k_D * norm(v) * v + obj.g];
        end

        function dfdx = dfdx(obj, t, x)
            v = x(4:6);
            v_norm = norm(v);
            dfdx = [zeros(3,3),  eye(3);
                    zeros(3,3),  (-obj.k_D * (v_norm * eye(3) + (v * v') / v_norm))];
        end

        function dfdu = dfdu(obj, t, x)
            dfdu = [zeros(3,3); 1/obj.m * eye(3)];
        end

        function add_impulse(obj)
            error('Impulsive control not implemented for Quadrotor')
        end

    end
end