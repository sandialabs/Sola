%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Clear workspace and add the SOLA source path.
clear;
close all;
clc;
rng(1342);                              % Random seed for reproducing results (optional).

%% Instantiate the objective and constraints.
T = 4 * pi;                             % Final simulation time.
n_t = 200;                              % Number of time steps.

radius = 1;
velocity = 1;
proportionality = radius^3 * velocity^2;

regularizer = 5;                        % Control regularization parameter.
r0 = 1.5;                               % Initial position (radius).
w0 = 2;                                 % Initial angular velocity.

objective = Tutorial_2B_Objective(T, n_t, radius, velocity, regularizer);
constraint = Tutorial_2B_Constraint(T, n_t, r0, w0, proportionality);

%% Solve the optimization problem with the Gauss-Newton Hessian approximation.
optimizer = Reduced_Space_Optimization(objective, constraint);
optimizer.Gauss_Newton_Hess = true;     % Use Algorithm 3, not Algorithm 2.

z0 = randn(objective.n_z, 1);           % Initial guess for the control.
[u, z] = optimizer.Optimize(z0);

%% Plot the results.
objective.Plot(u, z);
