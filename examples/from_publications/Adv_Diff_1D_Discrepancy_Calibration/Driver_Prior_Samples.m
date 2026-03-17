%% Set up
clear;
close all;
clc;
addpath(genpath('../../src'));
load Optimization_Results.mat;
rng(1232453);

working_path = pwd;
write_path = '~/Documents/dasco/papers/Model_Discrepancy_Hyperparameters/figures/';

alpha_u_pert = 100;
beta_u_pert = 100;
alpha_z_pert = 100;
beta_z_pert = 100;

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

u_hyperparam_interface = MD_u_Hyperparameter_Interface_Discrepancy_Calibration(x, data_centering);
u_prior_interface = MD_Numeric_Laplacian_u_Prior_Interface(con_hifi.S, con_hifi.M, data_interface, u_hyperparam_interface);

num_state_solves = 100;
z_hyperparam_interface = MD_z_Hyperparameter_Interface_Discrepancy_Calibration(num_state_solves, x, con_lofi);
z_prior_interface = MD_Numeric_Laplacian_z_Prior_Interface(con_hifi.S, con_hifi.M, data_interface, z_hyperparam_interface, u_prior_interface);

alpha_u = alpha_u_pert / 100 * u_hyperparam_interface.alpha_u;
u_prior_interface.Set_alpha_u(alpha_u);
beta_u = beta_u_pert / 100 * u_hyperparam_interface.beta_u;
u_prior_interface.Set_beta_u(beta_u);
u_vec = 0 * D(:, 1);
u_prior_interface.Compute_E_u_Inverse_GSVD(u_hyperparam_interface.gsvd_num_sing_vals, u_hyperparam_interface.gsvd_oversampling, u_hyperparam_interface.gsvd_num_subspace_iter, u_vec);

alpha_z = alpha_z_pert / 100 * z_hyperparam_interface.alpha_z;
z_prior_interface.Set_alpha_z(alpha_z);
beta_z = beta_z_pert / 100 * z_hyperparam_interface.beta_z;
z_prior_interface.Set_beta_z(beta_z);

%%
N = size(data_interface.Z, 2);
colors = distinguishable_colors(N);

num_prior_samples = 100;
md_prior_sampling = MD_Prior_Sampling(data_interface, u_prior_interface, z_prior_interface);

z = zeros(m, 2);
z(:, 1) = Z(:, 1) + .8 * x .* (1 - x);
z(:, 2) = Z(:, 1) .* (1 + .1 * cos(4 * pi * x));

figure;
hold on;
for k = 1:N
    plot(x, data_interface.Z(:, k), 'LineWidth', 3, 'color', colors(k, :));
end
plot(x, z(:, 1), 'LineWidth', 3, 'Color', 'magenta');
legend({'$z_1$', '$z_2$', '$z^{(1)}$'}, 'Interpreter', 'latex', 'Position', [0.7312    0.4914    0.1355    0.2204], 'FontSize', 24);
xlabel('Spatial Input');
ylabel('Source');
set(gca, 'fontsize', 24);
cd(write_path);
saveas(gca, 'z1_data', 'epsc');
cd(working_path);

figure;
hold on;
for k = 1:N
    plot(x, data_interface.Z(:, k), 'LineWidth', 3, 'color', colors(k, :));
end
plot(x, z(:, 2), 'LineWidth', 3, 'Color', 'cyan');
legend({'$z_1$', '$z_2$', '$z^{(2)}$'}, 'Interpreter', 'latex', 'Position', [0.7312    0.4914    0.1355    0.2204], 'FontSize', 24);
xlabel('Spatial Input');
ylabel('Source');
set(gca, 'fontsize', 24);
cd(write_path);
saveas(gca, 'z2_data', 'epsc');
cd(working_path);

figure;
hold on;
for k = 1:N
    plot(x, data_interface.D(:, k) + data_interface.data_shift, 'LineWidth', 3, 'color', colors(k, :));
end
legend({'$S(z_1)-\tilde{S}(z_1)$', '$S(z_2)-\tilde{S}(z_2)$'}, 'Interpreter', 'latex', 'Position', [0.2347    0.6604    0.2821    0.1682], 'FontSize', 24);
xlabel('Spatial Input');
ylabel('Discrepancy');
set(gca, 'fontsize', 24);
cd(write_path);
saveas(gca, 'delta_data', 'epsc');
cd(working_path);

[delta_samples, delta_zopt_samples] = md_prior_sampling.Prior_Discrepancy_Samples(z, num_prior_samples);

figure;
hold on;
plot(x, delta_zopt_samples, 'LineWidth', 3, 'color', [.9, .9, .9]);
plot(x, delta_zopt_samples(:, 1:6), 'LineWidth', 3);
xlabel('Spatial Input');
ylabel('Discrepancy');
ylim([-20, 15]);
set(gca, 'fontsize', 24);
cd(write_path);
saveas(gca, ['delta_0_samples_', num2str(alpha_u_pert), '_', num2str(beta_u_pert), '_', num2str(alpha_z_pert), '_', num2str(beta_z_pert)], 'epsc');
cd(working_path);

figure;
hold on;
for k = 1:100
    plot(x, delta_samples{k}(:, 1) - delta_zopt_samples(:, k), 'LineWidth', 3, 'Color', [.9, .9, .9]);
end
for k = 1:6
    plot(x, delta_samples{k}(:, 1) - delta_zopt_samples(:, k), 'LineWidth', 3);
end
xlabel('Spatial Input');
ylabel('Discrepancy');
ylim([-20, 15]);
set(gca, 'fontsize', 24);
cd(write_path);
saveas(gca, ['delta_z1_samples_', num2str(alpha_u_pert), '_', num2str(beta_u_pert), '_', num2str(alpha_z_pert), '_', num2str(beta_z_pert)], 'epsc');
cd(working_path);

figure;
hold on;
for k = 1:100
    plot(x, delta_samples{k}(:, 2) - delta_zopt_samples(:, k), 'LineWidth', 3, 'Color', [.9, .9, .9]);
end
for k = 1:6
    plot(x, delta_samples{k}(:, 2) - delta_zopt_samples(:, k), 'LineWidth', 3);
end
xlabel('Spatial Input');
ylabel('Discrepancy');
ylim([-20, 15]);
set(gca, 'fontsize', 24);
cd(write_path);
saveas(gca, ['delta_z2_samples_', num2str(alpha_u_pert), '_', num2str(beta_u_pert), '_', num2str(alpha_z_pert), '_', num2str(beta_z_pert)], 'epsc');
cd(working_path);
