clear;
close all;
clc;
addpath(genpath('../../src'));

%% Set up the optimization problem.
n_y = 200;
n_t = 51;
T = 1;
num_space_control_nodes = 10;
n_z = num_space_control_nodes * (n_t - 1);
obj = Adv_Diff_Gaussian_Source_Objective(n_y, n_z, T, n_t, num_space_control_nodes);
con = Adv_Diff_Gaussian_Source_Constraint(n_y, n_z, T, n_t, num_space_control_nodes);
opt = Reduced_Space_Optimization(obj, con);

%% Solve the optimization problem.
z0 = rand(n_z, 1);
[u, z] = opt.Optimize(z0);

%% Compare the optimal state to the target.
x = con.x;
t = con.t_mesh;

u_reshape = reshape(u, n_y, n_t);
figure;
for k = 1:n_t
    target = obj.Evaluate_Target(t(k), x);
    plot(x, u_reshape(:, k), '-', x, target, '--', 'LineWidth', 3);
    legend({'State', 'Target'});
    ylim([-1, 2]);
    pause(.05);
end
