% Clear Workspace and Add Interfaces to Path
clear;
close all;
addpath(genpath('../../src'));
% Retrieve Model Parameters (D, Z, diff/reg/react_coeff, m, u_lofi, z_hifi/lofi)
load Optimization_Results.mat;

% Set Default Font Axes and Line Width
set(0, "DefaultAxesFontSize", 20);
set(0, "DefaultLineLineWidth", 3);
set(0, "DefaultLineMarkerSize", 20);

% Set Hi-Fi and Lo-Fi Objectives and Constraints
obj = Diff_React_Objective(m, reg_coeff);
con_lofi = Diff_React_Constraint(m, diff_coeff, react_coeff);
opt_lofi = Reduced_Space_Optimization(obj, con_lofi);
con_hifi = Diff_React_HiFi_Constraint(con_lofi);
opt_hifi = Reduced_Space_Optimization(obj, con_hifi);
x = con_lofi.x;

%% HDSA interfaces
data_interface = MD_Data_Interface_Diff_React();
data_interface.Load_Data();

alpha_u = 2^2;
alpha_z = 1.e-10;
u_prior_interface = MD_Elliptic_u_Prior_Interface_Diff_React(alpha_u, opt_lofi);
z_prior_interface = MD_Elliptic_z_Prior_Interface_Diff_React(alpha_z, opt_lofi);

%% Hessian analysis
opt_prob_interface = MD_Opt_Prob_Interface_Sabl(opt_lofi, data_interface);
md_hessian_analysis = MD_Hessian_Analysis(opt_prob_interface, z_prior_interface);
num_evals = 4;
oversampling = 20;
md_hessian_analysis.Compute_Hessian_GEVP(data_interface.z_opt, num_evals, oversampling);

% Perform Offline OED Computations - USES data_interface
alpha_zd = 0.05;
beta_zd = 1.e-2;
oed_interface = MD_OED_Interface_Diff_React(data_interface, con_lofi, alpha_zd, beta_zd);
md_oed = MD_OED(opt_prob_interface, data_interface, u_prior_interface, z_prior_interface, md_hessian_analysis, oed_interface);
md_oed.Offline_Computation();

% Generate and Plot Random Design
N = 5;
Z = md_oed.Generate_Random_Design(N);
figure;
plot(x, Z, 'LineWidth', 3);
title('Random design');
set(gca, 'fontsize', 18);

md_oed.Offline_Computation();

% Generate and Plot Random Design from Subspace
Z = md_oed.Generate_Random_Design_from_Subspace(N);
figure;
plot(x, Z, 'LineWidth', 3);
title('Random design from subspace');
set(gca, 'fontsize', 18);

% Set Initial Guess for OED
beta_0 = randn(num_evals * (N - 1), 1);
alpha_d = 1.e-4;

% % Perform L-Curve Analysis to determine regularization parameter
% reg_coeffs = 10.^(-4:-1:-8)';
% [beta_L_curve, Z_L_curve, post_var, reg_val] = md_oed.L_Curve_Analysis(beta_0, alpha_d, reg_coeffs);
% figure;
% plot(post_var, reg_val, 'o', 'MarkerSize', 10);
% set(gca, 'fontsize', 18);
% % Choose the best result
% L_curve_index = 3;
% Z = Z_L_curve{L_curve_index};

reg_coeff = 1.e-6;
[betas, Z] = md_oed.Generate_Optimal_Design(beta_0, alpha_d, reg_coeff);

figure;
plot(x, Z, 'LineWidth', 3);
title('OED design');
set(gca, 'fontsize', 18);
