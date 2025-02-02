%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% minfuelparams.m returns a struct with parameters for the minfuel problem
% The parameters are used to define the problem in the SCPMinFuelSolver class
% Call this function to get the default parameters for the minfuel problem
% and modify them as needed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function p = minfuelparams()

    p.includeMass             = false; % include the mass variational equations for low-thrust
    p.CntrlType               = []; % 'impulsive' or 'LowThrust'
    p.DS                      = []; % DynamicalSystem class
    p.mesh_type               = 'fixed'; % 'fixed' or 'adaptive'
    p.x_init                  = []; % initial state
    p.x_fin                   = []; % final state
    p.tof                     = []; % time of flight
    p.u_min                   = 0; % minimum control magnitude
    p.u_max                   = []; % maximum control magnitude
    p.Nseg                    = []; % number of segments
    p.t_his                   = []; % time history for discretization
    p.s_min                   = 1E-6; % minimum time dilation for adaptive mesh
    p.s_max                   = []; % maximum time dilation for adaptive mesh

    % Collision avoidance
    p.n_obs                   = 0; % number of spherical obstacles
    p.R_obs                   = []; % radius of obstacles
    p.p_obs                   = []; % position of obstacles

    % Parameterized terminal state
    p.is_parameterized_x_init = false;
    p.is_parameterized_x_fin  = false;
    p.eval_x_init             = []; % callable function for evaluating the initial state with input tau_x_init
    p.eval_x_init_gradient    = []; % callable function for evaluating the gradient of the initial state with input tau_x_init
    p.eval_x_fin              = []; % callable function for evaluating the final state with input tau_x_fin
    p.eval_x_fin_gradient     = []; % callable function for evaluating the gradient of the final state with input tau_x_fin

end
