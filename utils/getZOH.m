function u = getZOH(t, t_his, u_his)
    
    % getZOH - get the value of the zero-order-hold control signal at time t
    %
    % Inputs:
    %    t - time at which to evaluate the control signal
    %    t_his - ZOH discretization time nodes (N+1 x 1)
    %    u_his - ZOH control signal values (n_u x N)

    % if t is outside the time span of the control history, return zero
    if t < t_his(1) || t > t_his(end)
        u = zeros(size(u_his, 1), 1);
        return;
    end
        
    % find the time segment right before t
    k = find(t_his <= t, 1, 'last');

    if isempty(k)
        error('t must be within the time span of the control history');
    end

    if k == length(t_his)
        u = u_his(:,end);
    else
        u = u_his(:,k);
    end
end
    