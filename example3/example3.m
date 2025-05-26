%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% example3: main script for demonstrating SCvx* for Earth-Mars transfer
% The problem is defined in EarthMars.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear

addpath(fullfile(fileparts(mfilename('fullpath')), '..')) % Path to parent folder with astro module
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'src')) % Path to SCvx* source code
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'utils')) % Path to utility functions

% load mice kernels
% Install the following files from the NAIF website into the kernels folder:
% --- de440.bsp, mar097.bsp, naif0012.tls
if ~exist(fullfile(fileparts(mfilename('fullpath')), '..', 'kernels', 'de440.bsp'), 'file') || ...
    ~exist(fullfile(fileparts(mfilename('fullpath')), '..', 'kernels', 'mar097.bsp'), 'file') || ...
    ~exist(fullfile(fileparts(mfilename('fullpath')), '..', 'kernels', 'naif0012.tls'), 'file')
    error('Please download the kernels from the NAIF website and place them in the kernels folder')
end
% Load the kernels
cspice_furnsh(fullfile(fileparts(mfilename('fullpath')), 'kernels', 'de440.bsp'))
cspice_furnsh(fullfile(fileparts(mfilename('fullpath')), 'kernels', 'mar097.bsp'))
cspice_furnsh(fullfile(fileparts(mfilename('fullpath')), 'kernels', 'naif0012.tls'))

% Define problem parameters
p = EarthMars();

% Setup initial guess
init_guess_struct.x = astro.conics.interp_MEE_from_cart(p.x_init, p.x_fin, p.Nseg+1, p.DS.mu, Nrev=0);

% Set up problem
prob = SCPMinFuelSolver(p, init_guess_struct);
prob.solve(save_bool = false);

%% Results
% Integrate the solution
prob.set_traj_with_stm();

%% Plot the trajectory
figure; hold on
plot3d(init_guess_struct.x, 'g--', 'DisplayName', 'Initial guess')
plot3d(prob.traj_with_stm.y, 'k', 'DisplayName', 'Solution')
labels3d('AU')
legend()

%% Plot the control profile
figure
stairsZOH(p.t_his*p.DS.SF.TU2DAYS, prob.sol.v*p.DS.SF.a*1E6)
xlabel('$t$ [days]')
ylabel('$\|u\|$ [mm/s$^2$]')
