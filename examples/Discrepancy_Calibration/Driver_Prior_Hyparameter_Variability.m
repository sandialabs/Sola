%% Set up
clear;
close all;
clc;
addpath(genpath('../../src'));
load Optimization_Results.mat;
rng(1232453);

working_path = pwd;
write_path = '~/Documents/dasco/papers/Model_Discrepancy_Hyperparameters/figures/';

con_hifi = Adv_Diff_Constraint(m, diff_coeff, vel_coeff, robin_coeff);
con_lofi = Diff_Constraint(con_hifi);
x = con_lofi.x;

K = 9;
M = 50;
N_range = 2:(K + 1);
alpha_u = zeros(K, M);
beta_u = zeros(K, M);
alpha_z = zeros(K, M);
beta_z = zeros(K, M);

E_z = (1.e-1) * con_hifi.S + con_hifi.M;
for N = 1:K
    num_samples = N_range(N);

    for j = 1:M
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

        alpha_u(N, j) = u_hyperparam_interface.alpha_u;
        beta_u(N, j) = u_hyperparam_interface.beta_u;
        alpha_z(N, j) = z_hyperparam_interface.alpha_z;
        beta_z(N, j) = z_hyperparam_interface.beta_z;
    end
end

disp('alpha_u');
disp(mean(alpha_u, 2)');
disp(std(alpha_u, [], 2)');
disp('beta_u');
disp(mean(beta_u, 2)');
disp(std(beta_u, [], 2)');
disp('alpha_z');
disp(mean(alpha_z, 2)');
disp(std(alpha_z, [], 2)');
disp('beta_z');
disp(mean(beta_z, 2)');
disp(std(beta_z, [], 2)');

colors = distinguishable_colors(K);

figure;
hold on;
q = quantile(alpha_u', [.0, 1]);
for k = 1:K
    plot([N_range(k), N_range(k)], [q(1, k), q(2, k)], 'Color', colors(k, :), 'LineWidth', 5);
end
plot(N_range, mean(alpha_u(:)) * .8 * ones(K, 1), '--', 'Color', 'cyan', 'LineWidth', 3);
plot(N_range, mean(alpha_u(:)) * 1.2 * ones(K, 1), '--', 'Color', 'cyan', 'LineWidth', 3);
xlabel('$N$', 'Interpreter', 'latex');
ylabel('$\alpha_u$', 'Interpreter', 'latex');
set(gca, 'fontsize', 24);
cd(write_path);
saveas(gca, 'alpha_u_variability', 'epsc');
cd(working_path);

figure;
hold on;
q = quantile(beta_u', [0, 1]);
for k = 1:K
    plot([N_range(k), N_range(k)], [q(1, k), q(2, k)], 'Color', colors(k, :), 'LineWidth', 5);
end
plot(N_range, mean(beta_u(:)) * .8 * ones(K, 1), '--', 'Color', 'cyan', 'LineWidth', 3);
plot(N_range, mean(beta_u(:)) * 1.2 * ones(K, 1), '--', 'Color', 'cyan', 'LineWidth', 3);
xlabel('$N$', 'Interpreter', 'latex');
ylabel('$\beta_u$', 'Interpreter', 'latex');
set(gca, 'fontsize', 24);
cd(write_path);
saveas(gca, 'beta_u_variability', 'epsc');
cd(working_path);

figure;
hold on;
q = quantile(alpha_z', [0, 1]);
for k = 1:K
    plot([N_range(k), N_range(k)], [q(1, k), q(2, k)], 'Color', colors(k, :), 'LineWidth', 5);
end
plot(N_range, mean(alpha_z(:)) * .8 * ones(K, 1), '--', 'Color', 'cyan', 'LineWidth', 3);
plot(N_range, mean(alpha_z(:)) * 1.2 * ones(K, 1), '--', 'Color', 'cyan', 'LineWidth', 3);
xlabel('$N$', 'Interpreter', 'latex');
ylabel('$\alpha_z$', 'Interpreter', 'latex');
set(gca, 'fontsize', 24);
cd(write_path);
saveas(gca, 'alpha_z_variability', 'epsc');
cd(working_path);

figure;
hold on;
q = quantile(beta_z', [0, 1]);
for k = 1:K
    plot([N_range(k), N_range(k)], [q(1, k), q(2, k)], 'Color', colors(k, :), 'LineWidth', 5);
end
plot(N_range, mean(beta_z(:)) * .8 * ones(K, 1), '--', 'Color', 'cyan', 'LineWidth', 3);
plot(N_range, mean(beta_z(:)) * 1.2 * ones(K, 1), '--', 'Color', 'cyan', 'LineWidth', 3);
xlabel('$N$', 'Interpreter', 'latex');
ylabel('$\beta_z$', 'Interpreter', 'latex');
ylim([0, .009]);
set(gca, 'fontsize', 24);
cd(write_path);
saveas(gca, 'beta_z_variability', 'epsc');
cd(working_path);
