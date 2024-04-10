%% Clear workspace and add the SABL source path.
clear;
close all;
clc;
run('~/Software/SABL/src/Set_Paths');
rng(1342);              % Random seed for reproducing exact results (optional).

%% Set the dimensions.
n_y = 4;
n_z = 2;
T = 2 * pi;
n_t = 100;
radius = 1;
velocity = 1;
proportionality = radius^3 * velocity^2;

%% Use automatic differentiation classes to solve the optimization problem.
obj_AD = Tutorial_2_Objective_AD(n_y, n_z, T, n_t, radius, velocity);
con_AD = Tutorial_2_Constraint_AD(n_y, n_z, T, n_t, proportionality);
opt_AD = Reduced_Space_Optimization(obj_AD, con_AD);

% evalc("obj_AD.AD_Initialization()");
% evalc("con_AD.AD_Initialization()");
obj_AD.AD_Initialization();
con_AD.AD_Initialization();

z0 = [0.9; 1.1];

[u, z] = opt_AD.Optimize(z0);

%% Post-process the results.
t = obj_AD.t_mesh;
ys = reshape(u, n_y, n_t);
r = ys(1, :);
theta = ys(3, :);
x = r .* cos(theta);
y = r .* sin(theta);

figure;
plot(t, x);
hold on;
plot(t, y);
title('x and y coordinates over time');
axis('equal');

figure;
plot(x, y);
title('position');
