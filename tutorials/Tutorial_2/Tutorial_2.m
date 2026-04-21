%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Clear workspace and add the SABL source path.
clear;
close all;
clc;
rng(18);                                % Random seed for reproducing results (optional).

%% Instantiate the objective and constraints.
T = 4 * pi;                             % Final simulation time.
n_t = 200;                              % Number of time steps.

radius = 1;
velocity = 1;
proportionality = radius^3 * velocity^2;

objective = Tutorial_2_Objective(T, n_t, radius, velocity);
constraint = Tutorial_2_Constraint(T, n_t, proportionality);

%% Solve the optimization problem.
optimizer = Reduced_Space_Optimization(objective, constraint);
z0 = 1 + randn(objective.n_z, 1) / 3;   % Initial guess for the control.
[u, z] = optimizer.Optimize(z0);

disp("Control:");
disp(z);

%% Plot the results.
objective.Plot(u);
