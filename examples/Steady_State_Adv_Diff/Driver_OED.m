%% Set up
clear;
close all;
clc;
addpath(genpath('../../src'));
load Optimization_Results.mat;

obj = Adv_Diff_Objective(m, reg_coeff);
con_hifi = Adv_Diff_Constraint(m, diff_coeff, vel_coeff, robin_coeff);
con_lofi = Diff_Constraint(con_hifi);
opt_hifi = Reduced_Space_Optimization(obj, con_hifi);
opt_lofi = Reduced_Space_Optimization(obj, con_lofi);
x = con_lofi.x;

%% HDSA interfaces
data_interface = MD_Data_Interface_Diff();
data_interface.Load_Data();

alpha_u = (1 / 2)^2;
alpha_z = (1 / 100)^2;
u_prior_interface = MD_Elliptic_u_Prior_Interface_Diff(alpha_u, opt_lofi);
z_prior_interface = MD_Elliptic_z_Prior_Interface_Diff(alpha_z, opt_lofi);

%% Hessian analysis
opt_prob_interface = MD_Opt_Prob_Interface_Sabl(opt_lofi, data_interface);
md_hessian_analysis = MD_Hessian_Analysis(opt_prob_interface, z_prior_interface);
num_evals = 20;
oversampling = 10;
md_hessian_analysis.Compute_Hessian_GEVP(data_interface.z_opt, num_evals, oversampling);

%% OED
oed_interface = MD_OED_Interface_Diff(data_interface, con_lofi);

md_oed = MD_OED(opt_prob_interface, data_interface, u_prior_interface, z_prior_interface, md_hessian_analysis, oed_interface);

N = 5;
Z = md_oed.Generate_Random_Design(N);
figure;
plot(x, Z, 'LineWidth', 3);
title('Random design');
set(gca, 'fontsize', 18);

md_oed.Offline_Computation();

Z = md_oed.Generate_Random_Design_from_Subspace(N);
figure;
plot(x, Z, 'LineWidth', 3);
title('Random design from subspace');
set(gca, 'fontsize', 18);

beta_0 = randn(num_evals * (N - 1), 1);
alpha_d = 1.e-2;
reg_coeffs = 10.^(-4:-1:-8)';
[beta_L_curve, Z_L_curve, post_var, reg_val] = md_oed.L_Curve_Analysis(beta_0, alpha_d, reg_coeffs);
figure;
plot(post_var, reg_val, 'o', 'MarkerSize', 10);
set(gca, 'fontsize', 18);

L_curve_index = 3;
Z = Z_L_curve{L_curve_index};

% reg_coeff = 1.e-6;
% [beta, Z] = md_oed.Generate_Optimal_Design(beta_0, alpha_d,reg_coeff);

figure;
plot(x, Z, 'LineWidth', 3);
title('OED design');
set(gca, 'fontsize', 18);
