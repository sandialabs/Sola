% Clear Workspace and Add Interfaces to Path
% clear;
close all;
% clc;
addpath(genpath('../../src'));
rng(0);

% Set Default Font Axes and Line Width
set(0, "DefaultAxesFontSize", 20);
set(0, "DefaultLineLineWidth", 3);
set(0, "DefaultLineMarkerSize", 20);

% Retrieve Model Parameters (D, Z, diff/reg/react_coeff, m, u_lofi, z_hifi/lofi)
load Optimization_Results.mat;

% Set Hi-Fi and Lo-Fi Objectives and Constraints
obj = Adv_Diff_Objective(m, reg_coeff);
con_hifi = Adv_Diff_Constraint(m, diff_coeff, vel_coeff, robin_coeff);
con_lofi = Diff_Constraint(con_hifi);
opt_hifi = Reduced_Space_Optimization(obj, con_hifi);
opt_lofi = Reduced_Space_Optimization(obj, con_lofi);
x = con_lofi.x;

% Note this doesn't contain access to Z/D yet.
data_interface = MD_Data_Interface_Diff(u_lofi, z_lofi);

% Generate Priors for u and z
alpha_u = (1 / 2)^2;
alpha_z = (1 / 100)^2;
alpha_d = 1.e-3;
u_prior_interface = MD_Elliptic_u_Prior_Interface_Diff(alpha_u, opt_lofi);
z_prior_interface = MD_Elliptic_z_Prior_Interface_Diff(alpha_z, opt_lofi);

% Calculate Relative OED Error (Lambda Function for now)
oed_z_error_fn = @(z) sqrt((z - z_hifi)' * z_prior_interface.Apply_M_z(z - z_hifi)) / sqrt(z_hifi' * z_prior_interface.Apply_M_z(z_hifi));

% %% Hessian analysis
opt_prob_interface = MD_Opt_Prob_Interface_Sabl(opt_lofi, data_interface);
md_hessian_analysis = MD_Hessian_Analysis(opt_prob_interface, z_prior_interface);
num_evals = 4;
oversampling = 20;
md_hessian_analysis.Compute_Hessian_GEVP(data_interface.z_opt, num_evals, oversampling);

% Perform Offline OED Computations - USES data_interface
oed_interface = MD_OED_Interface_Diff(data_interface, con_lofi);
md_oed = MD_OED(opt_prob_interface, data_interface, u_prior_interface, z_prior_interface, md_hessian_analysis, oed_interface);
md_oed.Offline_Computation();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kronecker product computation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Im = eye(m);
Mz = z_prior_interface.M;
% B1 = @(x) opt_prob_interface.Apply_Solution_Operator_z_Jacobian(opt_prob_interface.Apply_Misfit_Hessian([Im kron(Im, z_lofi' * Mz)] * x, u_lofi, z_lofi), z_lofi);
B1 = @(x) opt_prob_interface.Apply_Solution_Operator_z_Jacobian_Transpose(opt_prob_interface.Apply_Misfit_Hessian([Im kron(Im, z_lofi' * Mz)] * x, u_lofi, z_lofi), z_lofi);
% B2 = @(x) [zeros(m,m) reshape(Mz * reshape(x(m+1:end), m, m) * opt_prob_interface.Misfit_Gradient', m^2, 1)];
B2 = @(x) [zeros(m, m) kron(opt_prob_interface.Misfit_Gradient(u_lofi, z_lofi)', Mz)] * x;
B = @(x) B1(x) + B2(x);
PHinvB = @(x) md_hessian_analysis.Apply_Projected_RS_Hessian_Inverse(B(x));

% Test
% rng(0);
% tt_ex = randn(200*(200+1), 1);
% disp(size(PHinvB(tt_ex)))

theta_post_mean_tmp = md_update.Posterior_Theta_Mean_Temp();
z_update_this = z_lofi - PHinvB(theta_post_mean_tmp);
disp(norm(z_update_mean - (z_lofi - PHinvB_mean)) / norm(z_update_mean)); % About 5% undiagnosed error...
