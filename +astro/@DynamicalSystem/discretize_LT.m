function [Phi_A, Phi_B, Phi_c, trajs] = discretize_LT(obj, x_ref, u_ref, t_his)

    nx = obj.nx;
    nu = obj.nu;

    N = length(t_his) - 1;

    Phi_A = NaN(nx,nx,N);
    Phi_B = NaN(nx, nu, N);
    Phi_c = NaN(nx,1,N);
    trajs = cell(N,1);

    for k = 1:N
        [Phi_A(:,:,k), Phi_B(:,:,k), Phi_c(:,:,k), trajs{k}] = discretize_segment(obj, x_ref(:,k), u_ref(:,k), t_his(k:k+1));
    end

end

function [Phi_A, Phi_B, Phi_c, traj] = discretize_segment(obj, x_ref_k, u_ref_k, tspan)
    nx = obj.nx;
    nu = obj.nu;
    
    Y_k = [x_ref_k
        vec(eye(nx))
        zeros(nx * nu, 1)
        ];

    traj = obj.integrator(@(t, y) EOM_with_discrete_matrices(obj, t, y, u_ref_k), ...
                            tspan, ...
                            Y_k, ...
                            obj.odeopts);
    y_kp1 = traj.y(:,end);
    x_kp1 = y_kp1(1:nx);
    Phi_A = reshape(y_kp1(nx+1:nx+nx^2), [nx nx]);
    Phi_B = reshape(y_kp1(nx+nx^2+1:nx+nx^2+nx*nu), [nx nu]);
    Phi_c = x_kp1 - Phi_A * x_ref_k - Phi_B * u_ref_k;
end

function dYdt = EOM_with_discrete_matrices(obj, t, Y, u_k)

    x = Y(1:obj.nx);
    Phi_A = reshape(Y(obj.nx + 1:obj.nx + obj.nx ^ 2), [obj.nx, obj.nx]);
    Phi_B = reshape(Y(obj.nx + obj.nx ^ 2 + 1:obj.nx + obj.nx ^ 2 + obj.nx * obj.nu), [obj.nx, obj.nu]);

    dfdx = obj.dfdx(t, x);
    dfdu = obj.dfdu(t, x);

    dYdt = [
            obj.EOM_with_LT(t, x, u_k)
            vec(dfdx * Phi_A)
            vec(dfdx * Phi_B + dfdu)
            ];
end