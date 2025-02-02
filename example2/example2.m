%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Implementation of Example 2 from Oguri (2023)
% "Successive Convexification with Feasibility Guarantee via 
% Augmented  Lagrangian for Non-Convex Optimal Control Problems"
% The implementation of the problem is in ExampleClass2.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear;

addpath(fullfile(fileparts(mfilename('fullpath')), '..')) % Path to parent folder with astro module
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'src')) % Path to SCvx* source code
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'utils')) % Path to utility functions

prob = ExampleClass2();
% prob = ExampleClass2Optimizer();
SCvxParams = SCPParams();
tic; prob.solve(scp_params = SCvxParams, save_bool = false, clear_every_iter=true);toc

%% Plot the trajectory and control vectors
figure; axis equal;
hold on;

% Change the order of the states to plot the trajectory in the north-east-up frame
x_init_NWU = [prob.init_guess_struct.x(2,:); prob.init_guess_struct.x(3,:); prob.init_guess_struct.x(1,:)];
x_NWU = [prob.sol.x(2,:); prob.sol.x(3,:); prob.sol.x(1,:)];
u_NWU = [prob.sol.u(2,:); prob.sol.u(3,:); prob.sol.u(1,:)];

plot3d(x_init_NWU, 'r', DisplayName='Initial guess')
plot3d(x_NWU, 'k', DisplayName='Trajectory')
quiver3d(x_NWU(:,1:end-1), u_NWU, 'b')

up_min = min(x_NWU(3,:)) * 1.05;
up_max = max(x_NWU(3,:)) * 1.05;

[spX, spY, spZ] = cylinder(1);
sp1X = prob.p_obs(2,1) + spX * prob.R_obs(1);
sp1Y = prob.p_obs(3,1) + spY * prob.R_obs(1);
sp1Z = prob.p_obs(1,1) + spZ * (abs(up_min) + abs(up_max)) + up_min;
sp2X = prob.p_obs(2,2) + spX * prob.R_obs(2);
sp2Y = prob.p_obs(3,2) + spY * prob.R_obs(2);
sp2Z = prob.p_obs(1,2) + spZ * (abs(up_min) + abs(up_max)) + up_min;

p3 = surf(sp1X, sp1Y, sp1Z);
p4 = surf(sp2X, sp2Y, sp2Z);
p3.FaceAlpha = 0.1; p3.EdgeAlpha = 0.5;
p4.FaceAlpha = 0.1; p4.EdgeAlpha = 0.5;

xlabel('East, [m] ')
ylabel('North, [m]')
zlabel('Altitude, [m]')

%% Plot the control magnitude
figure;
hold on
stairsZOH(prob.t_his, prob.sol.Gamma, DisplayName='$\Gamma$')
stairsZOH(prob.t_his, vecnorm(prob.sol.u), Displayname='$\|u\|$')
yline([prob.T_min, prob.T_max], 'k--', HandleVisibility='off')
xlabel('time [s]')
ylabel('Thrust magnitude [N]')
legend()

%% Plot the SCvx* convergence history
prob.scp.plot_iter_history()

%% Plot the processing time of the algorithm
prob.scp.plot_time()