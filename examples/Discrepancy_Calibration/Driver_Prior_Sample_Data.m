%% Set up
clear;
close all;
clc;
addpath(genpath('../../src'));
load Optimization_Results.mat;
rng(1232453);

con_hifi = Adv_Diff_Constraint(m, diff_coeff, vel_coeff, robin_coeff);
con_lofi = Diff_Constraint(con_hifi);
x = con_lofi.x;

num_samples = 2;

E_z = (1.e-1) * con_hifi.S + con_hifi.M;

z = linsolve(E_z, sqrtm(con_hifi.M) * randn(m, num_samples - 1));
z_lofi_norm = sqrt(z_lofi' * con_hifi.M * z_lofi);
z_norm = sqrt(diag(z' * con_hifi.M * z));
for k = 1:(num_samples - 1)
    z(:, k) = z_lofi + .3 * z_lofi_norm * z(:, k) / z_norm(k);
end
Z = zeros(m, num_samples);
Z(:, 1) = z_lofi;
Z(:, 2:end) = z;
D = con_hifi.State_Solve(Z) - con_lofi.State_Solve(Z);

data_interface = MD_Data_Interface_Discrepancy_Calibration(z_lofi, u_lofi, Z, D);
data_centering = true;

u_hyperparams = MD_u_Hyperparameters_Discrepancy_Calibration(data_interface, x, data_centering);
u_prior_interface = MD_Analytic_Laplacian_u_Prior_Interface(con_hifi.M, u_hyperparams);

z_hyperparams = MD_z_Hyperparameters_Discrepancy_Calibration(data_interface, u_prior_interface, x);
z_prior_interface = MD_Numeric_Laplacian_z_Prior_Interface(con_hifi.S, con_hifi.M, z_hyperparams);

md_prior_sampling = MD_Prior_Sampling(data_interface, u_prior_interface, z_prior_interface);

num_prior_samples = 100;
num_perts = 10;
[delta_samples_z_opt, delta_samples_z_pert, z_pert] = md_prior_sampling.Prior_Discrepancy_Samples_for_Visualization(num_prior_samples, num_perts);

figure;
hold on;
plot(x, delta_samples_z_opt, 'LineWidth', 3, 'Color', [.9, .9, .9]);
plot(x, delta_samples_z_opt(:, 1:6), 'LineWidth', 3);
title('$\delta(z_{opt})$ samples', 'Interpreter', 'latex');

figure;
hold on;
plot(x, z_pert(:, 1), 'LineWidth', 3);
title('Leading perturbation $z_1$', 'Interpreter', 'latex');

figure;
hold on;
plot(x, delta_samples_z_pert{1}, 'LineWidth', 3, 'Color', [.9, .9, .9]);
plot(x, delta_samples_z_pert{1}(:, 1:6), 'LineWidth', 3);
title('$\delta(z_{1})-\delta(z_{opt})$ samples', 'Interpreter', 'latex');
