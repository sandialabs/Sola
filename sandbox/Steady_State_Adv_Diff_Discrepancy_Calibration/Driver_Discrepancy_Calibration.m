%% Set up
clear;
close all;
clc;
rng(231242343);
addpath(genpath('../../src'));
load('Optimization_Results_xi_1.mat');

con_hifi = Adv_Diff(m, vel_coeff, xi);
con_lofi = Diff(con_hifi);
x = con_lofi.x;

M = con_hifi.M;
S = con_hifi.S;

data_interface = MD_Data_Interface_Diff(xi);
data_interface.Load_Data();

u_hyperparam_interface = MD_u_Hyperparameter_Interface_Diff(x);
u_hyperparam_interface.Set_alpha_u(2.0);
u_hyperparam_interface.Set_beta_u(0.005);
u_prior_interface = MD_Numeric_Laplacian_u_Prior_Interface(S, M, data_interface, u_hyperparam_interface);

z_hyperparam_interface = MD_z_Hyperparameter_Interface_Diff(x);
z_hyperparam_interface.Set_alpha_z(2.0);
z_hyperparam_interface.Set_beta_z(0.008);
z_prior_interface = MD_Numeric_Laplacian_z_Prior_Interface(S, M, data_interface, z_hyperparam_interface, u_prior_interface);

md_post_sampling = MD_Posterior_Sampling(data_interface, u_prior_interface, z_prior_interface);
alpha_d = 1.e-5;
num_post_samples = 100;
md_post_sampling.Compute_Posterior_Data(alpha_d, num_post_samples);
Z_test = zeros(m, 3);
Z_test(:, 1:2) = Z;
Z_test(:, 3) = 0.5 * Z(:, 1) + 0.5 * Z(:, 2) + 0.7 * x;
[delta_mean_xi_1, delta_samples_xi_1] = md_post_sampling.Posterior_Discrepancy_Samples(Z_test);
dis_eval_xi_1 = con_hifi.State_Solve(Z_test) - con_lofi.State_Solve(Z_test);

%%
load('Optimization_Results_xi_0.99.mat');

con_hifi = Adv_Diff(m, vel_coeff, xi);
con_lofi = Diff(con_hifi);

data_interface = MD_Data_Interface_Diff(xi);
data_interface.Load_Data();

u_hyperparam_interface = MD_u_Hyperparameter_Interface_Diff(x);
u_hyperparam_interface.Set_alpha_u(2.0);
u_hyperparam_interface.Set_beta_u(0.005);
u_prior_interface = MD_Numeric_Laplacian_u_Prior_Interface(S, M, data_interface, u_hyperparam_interface);

z_hyperparam_interface = MD_z_Hyperparameter_Interface_Diff(x);
z_hyperparam_interface.Set_alpha_z(2.0);
z_hyperparam_interface.Set_beta_z(0.008);
z_prior_interface = MD_Numeric_Laplacian_z_Prior_Interface(S, M, data_interface, z_hyperparam_interface, u_prior_interface);

md_post_sampling = MD_Posterior_Sampling(data_interface, u_prior_interface, z_prior_interface);
alpha_d = 1.e-5;
md_post_sampling.Compute_Posterior_Data(alpha_d, num_post_samples);

[delta_mean_xi_2, delta_samples_xi_2] = md_post_sampling.Posterior_Discrepancy_Samples(Z_test);
dis_eval_xi_2 = con_hifi.State_Solve(Z_test) - con_lofi.State_Solve(Z_test);

%%
load('Optimization_Results_xi_0.8.mat');

con_hifi = Adv_Diff(m, vel_coeff, xi);
con_lofi = Diff(con_hifi);

data_interface = MD_Data_Interface_Diff(xi);
data_interface.Load_Data();

u_hyperparam_interface = MD_u_Hyperparameter_Interface_Diff(x);
u_hyperparam_interface.Set_alpha_u(2.0);
u_hyperparam_interface.Set_beta_u(0.005);
u_prior_interface = MD_Numeric_Laplacian_u_Prior_Interface(S, M, data_interface, u_hyperparam_interface);

z_hyperparam_interface = MD_z_Hyperparameter_Interface_Diff(x);
z_hyperparam_interface.Set_alpha_z(2.0);
z_hyperparam_interface.Set_beta_z(0.008);
z_prior_interface = MD_Numeric_Laplacian_z_Prior_Interface(S, M, data_interface, z_hyperparam_interface, u_prior_interface);

md_post_sampling = MD_Posterior_Sampling(data_interface, u_prior_interface, z_prior_interface);
alpha_d = 1.e-5;
md_post_sampling.Compute_Posterior_Data(alpha_d, num_post_samples);

[delta_mean_xi_3, delta_samples_xi_3] = md_post_sampling.Posterior_Discrepancy_Samples(Z_test);
dis_eval_xi_3 = con_hifi.State_Solve(Z_test) - con_lofi.State_Solve(Z_test);

%%
r = max([max(max(abs(delta_samples_xi_1{3} - delta_samples_xi_3{3}))), max(max(abs(delta_samples_xi_1{3} - delta_samples_xi_2{3})))]) * 1.02;

figure;
hold on;
plot(x, dis_eval_xi_1(:, 1) - dis_eval_xi_2(:, 1), "Color", 'magenta', 'LineWidth', 3);
plot(x, delta_samples_xi_1{1} - delta_samples_xi_2{1}, "Color", [.9, .9, .9], 'LineWidth', 3);
plot(x, dis_eval_xi_1(:, 1) - dis_eval_xi_2(:, 1), "Color", 'magenta', 'LineWidth', 3);
legend({'$(S(z_1,\xi_1)-\tilde{S}(z_1,\xi_1))-(S(z_1,\xi_2)-\tilde{S}(z_1,\xi_2))$', '$\delta(z_1,\xi_1)-\delta(z_1,\xi_2)$'}, 'Interpreter', 'latex');
ylim([-r, r]);
set(gca, 'fontsize', 20);

figure;
hold on;
plot(x, dis_eval_xi_1(:, 3) - dis_eval_xi_2(:, 3), "Color", 'magenta', 'LineWidth', 3);
plot(x, delta_samples_xi_1{3} - delta_samples_xi_2{3}, "Color", [.9, .9, .9], 'LineWidth', 3);
plot(x, dis_eval_xi_1(:, 3) - dis_eval_xi_2(:, 3), "Color", 'magenta', 'LineWidth', 3);
legend({'$(S(z_3,\xi_1)-\tilde{S}(z_3,\xi_1))-(S(z_3,\xi_2)-\tilde{S}(z_3,\xi_2))$', '$\delta(z_3,\xi_1)-\delta(z_3,\xi_2)$'}, 'Interpreter', 'latex');
ylim([-r, r]);
set(gca, 'fontsize', 20);

figure;
hold on;
plot(x, dis_eval_xi_1(:, 1) - dis_eval_xi_3(:, 1), "Color", 'magenta', 'LineWidth', 3);
plot(x, delta_samples_xi_1{1} - delta_samples_xi_3{1}, "Color", [.9, .9, .9], 'LineWidth', 3);
plot(x, dis_eval_xi_1(:, 1) - dis_eval_xi_3(:, 1), "Color", 'magenta', 'LineWidth', 3);
legend({'$(S(z_1,\xi_1)-\tilde{S}(z_1,\xi_1))-(S(z_1,\xi_3)-\tilde{S}(z_1,\xi_3))$', '$\delta(z_1,\xi_1)-\delta(z_1,\xi_3)$'}, 'Interpreter', 'latex');
ylim([-r, r]);
set(gca, 'fontsize', 20);

figure;
hold on;
plot(x, dis_eval_xi_1(:, 3) - dis_eval_xi_3(:, 3), "Color", 'magenta', 'LineWidth', 3);
plot(x, delta_samples_xi_1{3} - delta_samples_xi_3{3}, "Color", [.9, .9, .9], 'LineWidth', 3);
plot(x, dis_eval_xi_1(:, 3) - dis_eval_xi_3(:, 3), "Color", 'magenta', 'LineWidth', 3);
legend({'$(S(z_3,\xi_1)-\tilde{S}(z_3,\xi_1))-(S(z_3,\xi_3)-\tilde{S}(z_3,\xi_3))$', '$\delta(z_3,\xi_1)-\delta(z_3,\xi_3)$'}, 'Interpreter', 'latex');
ylim([-r, r]);
set(gca, 'fontsize', 20);
