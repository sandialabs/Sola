%% Clear workspace and add the SABL source path.
clear;
close all;
clc;
run('~/Software/SABL/src/Set_Paths');
rng(1342);              % Random seed for reproducing exact results (optional).

%% Set dimensions and problem parameters.
n_y = 4;
n_q = 2;
T = 4 * pi;
n_t = 100;

radius = 1;
velocity = 1;
proportionality = radius^3 * velocity^2;
r0 = 1.5;
w0 = 3;

%% Solve the optimization problem.
objective = Tutorial_2B_Objective(T, n_t, radius, velocity);
constraint = Tutorial_2B_Constraint(T, n_t, r0, w0, proportionality);
optimization = Reduced_Space_Optimization(objective, constraint);

z0 = randn(objective.n_z, 1);   % Initial guess for the control.

[u, z] = optimization.Optimize(z0);

%% Plot the results.
objective.Plot(u, z);
