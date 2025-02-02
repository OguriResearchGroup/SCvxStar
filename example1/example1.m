%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implementation of Example 1 from Oguri (2023)
% "Successive Convexification with Feasibility Guarantee via 
% Augmented  Lagrangian for Non-Convex Optimal Control Problems"
% The implementation of the problem is in ExampleClass.m
% ExampleClass.m should be used as a template for creating new SCP problems.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; 
clear;

addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'src')) % Path to SCvx* source code
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'utils')) % Path to utility functions

SCvxParams = SCPParams();
% Set parameters to be the same as in the paper
SCvxParams.w_max       = 1E8;
SCvxParams.tol_opt     = 1E-5;
SCvxParams.tol_feas    = 1E-5;
SCvxParams.rho0        = 0;
SCvxParams.rho1        = 0.25;
SCvxParams.rho2        = 0.7;
SCvxParams.alpha1      = 2;
SCvxParams.alpha2      = 3;
SCvxParams.beta        = 2;
SCvxParams.gamma       = 0.9;
SCvxParams.r_init      = 0.1;
SCvxParams.r_min       = 1E-10;
SCvxParams.r_max       = 10;
SCvxParams.superlinear = false;

% Solve with different initial penalty values
% for w_init = [1E-1 1 1E1 1E2 1E3 1E4 1E5]

w_init = 1E3;
disp("Solving with w_init = " + w_init)
SCvxParams.w_init      = w_init;

% Define the objective function and constraints in a separate class
% Use yalmip's optimize (overhead at each iteration for problem construction)
prob = ExampleClass1();
% Use yalmip's optimizer (precompiles the problem at the beginning)
% prob = ExampleClass1Optimizer();
prob.solve(scp_params=SCvxParams, verbose=true, save_bool=false);

% Print the results
iters = prob.scp.report.iters;
fprintf('w_init = %e, iters = %d\n', w_init, iters)

% end
%% Plot feasible region and solution
figure
prob.plot_contour()
plot(prob.sol.z(1), prob.sol.z(2), 'b*')

%% Plot the convergence of the algorithm
prob.scp.plot_iter_history()

%% Plot the processing time of the algorithm
prob.scp.plot_time()