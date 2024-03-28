clear;
close all;
clc;
addpath(genpath('../../src'));

%% Set up the optimization problem.
n_y = 100;
n_t = 51;
T = 1;
num_space_control_nodes = 10;
reg_coeff = 1.e-6;
diff_coeff = 1;
vel_coeff_lofi = 1;
vel_coeff_hifi = 5;

n_z = num_space_control_nodes * (n_t - 1);
obj = Adv_Diff_Gaussian_Source_Objective(n_y, n_z, T, n_t, num_space_control_nodes, reg_coeff);

con_lofi = Adv_Diff_Gaussian_Source_Constraint(n_y, n_z, T, n_t, num_space_control_nodes, diff_coeff, vel_coeff_lofi);
opt_lofi = Reduced_Space_Optimization(obj, con_lofi);

con_hifi = Adv_Diff_Gaussian_Source_Constraint(n_y, n_z, T, n_t, num_space_control_nodes, diff_coeff, vel_coeff_hifi);
opt_hifi = Reduced_Space_Optimization(obj, con_hifi);

%% Solve the optimization problem.
z0 = rand(n_z, 1);
[u_lofi, z_lofi] = opt_lofi.Optimize(z0);
[u_hifi, z_hifi] = opt_hifi.Optimize(z_lofi);

%% Compare the optimal state to the target.
x = con_lofi.x;
t = con_lofi.t_mesh;
u_lofi_eval = con_hifi.State_Solve(z_lofi);

lb_lofi = min(u_lofi_eval) - abs(min(u_lofi_eval)) * .05;
ub_lofi = max(u_lofi_eval) + abs(max(u_lofi_eval)) * .05;
u_lofi_reshape = reshape(u_lofi_eval, n_y, n_t);

lb_hifi = min(u_hifi) - abs(min(u_hifi)) * .05;
ub_hifi = max(u_hifi) + abs(max(u_hifi)) * .05;
u_hifi_reshape = reshape(u_hifi, n_y, n_t);

lb = min(lb_lofi, lb_hifi);
ub = max(ub_lofi, ub_hifi);

figure;
for k = 1:n_t
    target = obj.Evaluate_Target(t(k), x);
    plot(x, u_hifi_reshape(:, k), '-', x, u_lofi_reshape(:, k), '-', x, target, '--', 'LineWidth', 3);
    legend({'Hifi State', 'Lofi State', 'Target'});
    ylim([lb, ub]);
    pause(.05);
end

Z = z_lofi;
D = u_lofi_eval - u_lofi;

save('Optimization_Results.mat', 'n_y', 'n_t', 'T', 'num_space_control_nodes', 'diff_coeff', 'vel_coeff_lofi', 'vel_coeff_hifi', 'reg_coeff', 'z_lofi', 'u_lofi', 'z_hifi', 'u_hifi', 'Z', 'D');
