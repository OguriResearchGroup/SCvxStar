function [A_, B_, c_, d_] = discretize_impulsive_adaptive(obj, x_ref_postDV, s_ref, tau_his)

    nx = obj.nx;
    N = length(tau_his) - 1;
    t_his = [0; cumsum(s_ref)/N];

    A_ = zeros(nx,nx,N);
    c_ = zeros(nx,1,N);
    d_ = zeros(nx,1,N);

    parfor k = 1:N
        [A_(:,:,k), c_(:,:,k), d_(:,:,k)] ...
            = discretize_segment(obj, x_ref_postDV(:,k), s_ref(k), tau_his(k:k+1), tau_his, t_his);
    end

    B_ = repmat(obj.B_impulse, [1 1 N+1]);

end

function [A, c, d] = discretize_segment(obj, x_ref_k_postDV, s_ref_k, tau_span, tau_his, t_his)
    nx = obj.nx;
    
    Y = [x_ref_k_postDV
        vec(eye(nx))
        zeros(nx, 1)
        zeros(nx, 1)
        ];

    [~, y] = obj.integrator(@(tau, y) EOM_with_discrete_matrices_adaptive(obj, tau, y, s_ref_k, tau_his, t_his), ...
                            tau_span, ...
                            Y, ...
                            obj.odeopts);
    y_f = y(end,:)';
    [A, c, d] = unpack_matrices(y_f, nx);
end

function dYdt = EOM_with_discrete_matrices_adaptive(obj, tau, Y, s_k, tau_his, t_his)

    nx = obj.nx;
    
    % Convert from nondimensional time to dimensional time
    t = tau2t(tau, tau_his, t_his);

    x = Y(1:nx);
    STM = reshape(Y(nx + 1:nx + nx ^ 2), [nx, nx]);

    A = s_k * obj.dfdx(t, x);
    c = - A * x;
    d = obj.EOM(t, x);

    dYdt = [
            s_k * obj.EOM(t, x)
            vec(A * STM)
            STM \ c
            STM \ d
            ];
end
    
function [A, c, d] = unpack_matrices(y, nx)
    A = reshape(y(nx+1:nx+nx^2), [nx nx]);
    c = A * y(nx+nx^2+1:2*nx+nx^2);
    d = A * y(2 * nx + nx ^ 2 + 1:3 * nx + nx ^ 2);
end

function t = tau2t(tau, tau_his, t_his)
    t = interp1(tau_his, t_his, tau);
end