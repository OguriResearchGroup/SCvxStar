function [Phi_A, Phi_B, Phi_c, trajs] = discretize_impulsive(obj, x_ref_postDV, t_his)

    nx = obj.nx;
    N = length(t_his) - 1;

    Phi_A = NaN(nx,nx,N);
    Phi_c = NaN(nx,1,N);
    trajs = cell(N,1);

    parfor k = 1:N
        [Phi_A(:,:,k), Phi_c(:,:,k), trajs{k}] = discretize_segment(obj, x_ref_postDV(:,k), t_his(k:k+1));
    end

    Phi_B = repmat(obj.B_impulse, [1 1 N+1]);
end

function [Phi_A, Phi_c, traj] = discretize_segment(obj, x_ref_k_postDV, tspan)
    nx = obj.nx;
    
    Y_k = [x_ref_k_postDV
        vec(eye(nx))
        zeros(nx, 1)
        ];

    traj = obj.integrator(@(t, y) EOM_with_discrete_matrices(obj, t, y), ...
                            tspan, ...
                            Y_k, ...
                            obj.odeopts);
    y_f = traj.y(:,end);
    x_kp1 = y_f(1:nx);
    Phi_A = reshape(y_f(nx+1:nx+nx^2), [nx nx]);
    Phi_c = x_kp1 - Phi_A * x_ref_k_postDV;
end

function dYdt = EOM_with_discrete_matrices(obj, t, Y)

    x = Y(1:obj.nx);
    Phi_A = reshape(Y(obj.nx + 1:obj.nx + obj.nx ^ 2), [obj.nx, obj.nx]);

    dfdx = obj.dfdx(t, x);

    dYdt = [
            obj.EOM(t, x)
            vec(dfdx * Phi_A)
            ];
end