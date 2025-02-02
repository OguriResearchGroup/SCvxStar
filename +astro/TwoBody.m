classdef TwoBody < astro.DynamicalSystem
    properties
        mu % Gravitational parameter
    end

    properties (Constant)
        is_autonomous = true;
        is_control_affine = true;
        nx = 6;
        nu = 3;
        B_impulse = [zeros(3); eye(3)]; % Control matrix for impulsive control
    end

    methods
        function obj = TwoBody(mu, SF)
            obj@astro.DynamicalSystem(SF);
            obj.mu = mu;
        end

        function dxdt = EOM(obj, t, x)
            r = x(1:3);
            dxdt = [x(4:6); -obj.mu/norm(r)^3*r];
        end

        function A = dfdx(obj, t, x)
            r = x(1:3);
            A = [zeros(3) eye(3); obj.mu/norm(r)^3 * (3/norm(r)^2*(r*r')-eye(3))  zeros(3)];
        end

        function dfdu = dfdu(obj, t, x)
            dfdu = [zeros(3); eye(3)];
        end

        function x_plus = add_impulse(obj, x_minus, u)
            switch size(x_minus, 1)
                case 6
                    x_plus = x_minus + [zeros(3,1); u];
                case 42
                    x_plus = x_minus + [zeros(3,1); u; zeros(36, 1)];
                otherwise
                    error('Invalid state dimension. Must be 6 or 42.');
            end
        end

    end

end
