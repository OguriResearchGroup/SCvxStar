function [A_, B_, c_, d_] = discretize_LT_adaptive(obj, x_ref, u_ref, s_ref, tau_his)

    nx = obj.nx;
    nu = obj.nu;

    N = length(tau_his) - 1;
    t_his = [0; cumsum(s_ref)/N];

    A_ = zeros(nx,nx,N);
    B_ = zeros(nx, nu, N);
    c_ = zeros(nx,1,N);
    d_ = zeros(nx,1,N);

    parfor k = 1:N
        [A_(:,:,k), B_(:,:,k), c_(:,:,k), d_(:,:,k)] ...
            = discretize_segment(obj, x_ref(:,k), u_ref(:,k), s_ref(k), tau_his(k:k+1), tau_his, t_his);
    end

end

function [A, B, c, d] = discretize_segment(obj, x_ref_k, u_ref_k, s_ref_k, tau_span, tau_his, t_his)
    nx = obj.nx;
    nu = obj.nu;
    
    Y = [x_ref_k
        vec(eye(nx))
        zeros(nx * nu, 1)
        zeros(nx, 1)
        zeros(nx, 1)
        ];

    [~, y] = obj.integrator(@(tau, y) EOM_with_LT_with_discrete_matrices_adaptive(obj, tau, y, u_ref_k, s_ref_k, tau_his, t_his), ...
                            tau_span, ...
                            Y, ...
                            obj.odeopts);
    y_f = y(end,:)';
    [A, B, c, d] = unpack_matrices(y_f, nx, nu);
end

function dYdt = EOM_with_LT_with_discrete_matrices_adaptive(obj, tau, Y, u_k, s_k, tau_his, t_his)

    nx = obj.nx;
    
    % Convert from nondimensional time to dimensional time
    t = tau2t(tau, tau_his, t_his);

    x = Y(1:nx);
    STM = reshape(Y(nx + 1:nx + nx ^ 2), [nx, nx]);

    A = s_k * obj.dfdx(t, x);
    B = s_k * obj.dfdu(t, x);
    c = - A * x - B * u_k;
    d = obj.EOM_with_LT(t, x, u_k);

    dYdt = [
            s_k * obj.EOM_with_LT(t, x, u_k)
            vec(A * STM)
            vec(STM \ B)
            STM \ c
            STM \ d
            ];
end
    
function [A, B, c, d] = unpack_matrices(y, nx, nu)
    A = reshape(y(nx+1:nx+nx^2), [nx nx]);
    B = A * reshape(y(nx+nx^2+1:nx+nx^2+nx*nu), [nx nu]);
    c = A * y(nx+nx^2+nx*nu+1:2*nx+nx^2+nx*nu);
    d = A * y(2 * nx + nx ^ 2 + nx * nu + 1:3 * nx + nx ^ 2 + nx * nu);
end

function t = tau2t(tau, tau_his, t_his)
    t = interp1(tau_his, t_his, tau);
end