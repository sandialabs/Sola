%% Set up
clear;
close all;
clc;
rng(1242343)
addpath(genpath('../../src'));
load('Optimization_Results_xi_1.mat');

obj = Adv_Diff_Objective(m, reg_coeff);
con_hifi = Adv_Diff(m, diff_coeff, vel_coeff, xi);
con_lofi = Diff(con_hifi);
opt_hifi = Reduced_Space_Optimization(obj, con_hifi);
opt_lofi = Reduced_Space_Optimization(obj, con_lofi);
x = con_lofi.x;

M = con_hifi.M;
S = con_hifi.S;

data_interface = MD_Data_Interface_Diff(xi);
data_interface.Load_Data();

u_hyperparam_interface = MD_u_Hyperparameter_Interface_Diff(x);
u_prior_interface = MD_Numeric_Laplacian_u_Prior_Interface(S, M, data_interface, u_hyperparam_interface);


z_hyperparam_interface = MD_z_Hyperparameter_Interface_Diff(x);
z_prior_interface = MD_Numeric_Laplacian_z_Prior_Interface(S, M, data_interface, z_hyperparam_interface, u_prior_interface);

md_post_sampling = MD_Posterior_Sampling(data_interface, u_prior_interface, z_prior_interface);
alpha_d = u_hyperparam_interface.alpha_d;
num_post_samples = 100;
md_post_sampling.Compute_Posterior_Data(alpha_d, num_post_samples);
Z_test = zeros(m, 3);
Z_test(:, 1:2) = Z;
Z_test(:, 3) = 0.5 * Z(:,1) + 0.5 * Z(:,2) + 0.7 * x;
[delta_mean_xi_1, delta_samples_xi_1] = md_post_sampling.Posterior_Discrepancy_Samples(Z_test);
dis_eval_xi_1 = con_hifi.State_Solve(Z_test) - con_lofi.State_Solve(Z_test);


%%
load('Optimization_Results_xi_0.99.mat');

obj = Adv_Diff_Objective(m, reg_coeff);
con_hifi = Adv_Diff(m, diff_coeff, vel_coeff, xi);
con_lofi = Diff(con_hifi);
opt_hifi = Reduced_Space_Optimization(obj, con_hifi);
opt_lofi = Reduced_Space_Optimization(obj, con_lofi);

data_interface = MD_Data_Interface_Diff(xi);
data_interface.Load_Data();

u_hyperparam_interface = MD_u_Hyperparameter_Interface_Diff(x);
u_prior_interface = MD_Numeric_Laplacian_u_Prior_Interface(S, M, data_interface, u_hyperparam_interface);

z_hyperparam_interface = MD_z_Hyperparameter_Interface_Diff(x);
z_prior_interface = MD_Numeric_Laplacian_z_Prior_Interface(S, M, data_interface, z_hyperparam_interface, u_prior_interface);

md_post_sampling = MD_Posterior_Sampling(data_interface, u_prior_interface, z_prior_interface);
alpha_d = u_hyperparam_interface.alpha_d;
md_post_sampling.Compute_Posterior_Data(alpha_d, num_post_samples);

[delta_mean_xi_99, delta_samples_xi_99] = md_post_sampling.Posterior_Discrepancy_Samples(Z_test);
dis_eval_xi_99 = con_hifi.State_Solve(Z_test) - con_lofi.State_Solve(Z_test);

%%
load('Optimization_Results_xi_0.01.mat');

obj = Adv_Diff_Objective(m, reg_coeff);
con_hifi = Adv_Diff(m, diff_coeff, vel_coeff, xi);
con_lofi = Diff(con_hifi);
opt_hifi = Reduced_Space_Optimization(obj, con_hifi);
opt_lofi = Reduced_Space_Optimization(obj, con_lofi);

data_interface = MD_Data_Interface_Diff(xi);
data_interface.Load_Data();

u_hyperparam_interface = MD_u_Hyperparameter_Interface_Diff(x);
u_prior_interface = MD_Numeric_Laplacian_u_Prior_Interface(S, M, data_interface, u_hyperparam_interface);

z_hyperparam_interface = MD_z_Hyperparameter_Interface_Diff(x);
z_prior_interface = MD_Numeric_Laplacian_z_Prior_Interface(S, M, data_interface, z_hyperparam_interface, u_prior_interface);

md_post_sampling = MD_Posterior_Sampling(data_interface, u_prior_interface, z_prior_interface);
alpha_d = u_hyperparam_interface.alpha_d;
md_post_sampling.Compute_Posterior_Data(alpha_d, num_post_samples);

[delta_mean_xi_01, delta_samples_xi_01] = md_post_sampling.Posterior_Discrepancy_Samples(Z_test);
dis_eval_xi_01 = con_hifi.State_Solve(Z_test) - con_lofi.State_Solve(Z_test);

%%
r = max( [ max(max(abs(delta_samples_xi_1{3}-delta_samples_xi_01{3}))) , max(max(abs(delta_samples_xi_1{3}-delta_samples_xi_99{3}))) ]) * 1.02;

figure,
hold on
plot(x,delta_samples_xi_1{1}-delta_samples_xi_99{1},"Color",[.9,.9,.9],'LineWidth',3)
plot(x,dis_eval_xi_1(:,1)-dis_eval_xi_99(:,1),"Color",'magenta','LineWidth',3)
title('$\delta(z_1,\xi_1)-\delta(z_1,\xi_2)$','Interpreter','latex')
ylim([-r,r])
set(gca, 'fontsize', 20);

figure,
hold on
plot(x,delta_samples_xi_1{3}-delta_samples_xi_99{3},"Color",[.9,.9,.9],'LineWidth',3)
plot(x,dis_eval_xi_1(:,3)-dis_eval_xi_99(:,3),"Color",'magenta','LineWidth',3)
title('$\delta(z_3,\xi_1)-\delta(z_3,\xi_2)$','Interpreter','latex')
ylim([-r,r])
set(gca, 'fontsize', 20);

figure,
hold on
plot(x,delta_samples_xi_1{1}-delta_samples_xi_01{1},"Color",[.9,.9,.9],'LineWidth',3)
plot(x,dis_eval_xi_1(:,1)-dis_eval_xi_01(:,1),"Color",'magenta','LineWidth',3)
title('$\delta(z_1,\xi_1)-\delta(z_1,\xi_3)$','Interpreter','latex')
ylim([-r,r])
set(gca, 'fontsize', 20);

figure,
hold on
plot(x,delta_samples_xi_1{3}-delta_samples_xi_01{3},"Color",[.9,.9,.9],'LineWidth',3)
plot(x,dis_eval_xi_1(:,3)-dis_eval_xi_01(:,3),"Color",'magenta','LineWidth',3)
title('$\delta(z_3,\xi_1)-\delta(z_3,\xi_3)$','Interpreter','latex')
ylim([-r,r])
set(gca, 'fontsize', 20);