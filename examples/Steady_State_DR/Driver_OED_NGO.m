% Clear Workspace and Add Interfaces to Path
clear;
close all;
% clc;
addpath(genpath('../../src'));

% Set Default Font Axes and Line Width
set(0, "DefaultAxesFontSize", 20);
set(0, "DefaultLineLineWidth", 3);
set(0, "DefaultLineMarkerSize", 20);

% Retrieve Model Parameters (D, Z, diff/reg/react_coeff, m, u_lofi, z_hifi/lofi)
load Optimization_Results.mat;

% Set Hi-Fi and Lo-Fi Objectives and Constraints
obj = Diff_React_Objective(m, reg_coeff);
con_lofi = Diff_React_Constraint(m, diff_coeff, react_coeff);
opt_lofi = Reduced_Space_Optimization(obj, con_lofi);
con_hifi = Diff_React_HiFi_Constraint(con_lofi);
opt_hifi = Reduced_Space_Optimization(obj, con_hifi);
x = con_lofi.x;

% TODO: Allow for better direct loading of data from files (using load("", "").(""))
% Note this doesn't contain access to Z/D yet.
data_interface = MD_Data_Interface_Diff_React(u_lofi, z_lofi);

% Generate Priors for u and z
alpha_u = 2;
alpha_z = 1e-8;
alpha_d = (1.e-2)^2 * alpha_u;
u_prior_interface = MD_Elliptic_u_Prior_Interface_Diff_React(alpha_u, opt_lofi);
z_prior_interface = MD_Elliptic_z_Prior_Interface_Diff_React(alpha_z, opt_lofi);

% Calculate Relative OED Error (Lambda Function for now)
oed_z_error_fn = @(z) sqrt((z - z_hifi)' * z_prior_interface.Apply_M_z(z - z_hifi)) / sqrt(z_hifi' * z_prior_interface.Apply_M_z(z_hifi));

% %% Hessian analysis
opt_prob_interface = MD_Opt_Prob_Interface_Sabl(opt_lofi, data_interface);
md_hessian_analysis = MD_Hessian_Analysis(opt_prob_interface, z_prior_interface);
num_evals = 4;
oversampling = 20;
md_hessian_analysis.Compute_Hessian_GEVP(data_interface.z_opt, num_evals, oversampling);

% Perform Offline OED Computations - USES data_interface
alpha_zd = 1.e-2;
beta_zd = 1.e-2;
oed_interface = MD_OED_Interface_Diff_React(data_interface, con_lofi, alpha_zd, beta_zd);
md_oed = MD_OED_NGO(opt_prob_interface, data_interface, u_prior_interface, z_prior_interface, md_hessian_analysis, oed_interface);
md_oed.Offline_Computation();

% Set Parameters for OED
N = 5;
rng(0);
beta_0 = randn(num_evals * (N - 1), 1);
reg_coeff = 1.e-6;
% [betas, Z] = md_oed.Generate_Optimal_Design(beta_0, alpha_d, reg_coeff);

% Finite difference check...
% In = eye(length(beta_0));
% h = 1.e-6;
% [val1, grad] = md_oed.Evaluate_Posterior_Cov_Trace(beta_0, alpha_d);
% disp(grad'*In(:,1))
% [val2, ~] = md_oed.Evaluate_Posterior_Cov_Trace(beta_0 + h*In(:, 1), alpha_d);
% disp(1/h * (val2-val1))

% Generate Design (Generate_Random_Design(N), Generate_Random_Design_from_Subspace(N), Generate_Optimal_Design(...))
% Z = md_oed.Generate_Random_Design(N);
D = Evaluate_Discrepancy(con_hifi, con_lofi, Z);
data_interface.Set_Z_and_D(Z, D);

% % Sample from Posterior (i.e., solve problem) - USES DATA INTERFACE
num_post_samples = 1;
md_post_sampling = MD_Posterior_Sampling(data_interface, u_prior_interface, z_prior_interface);
md_post_sampling.Compute_Posterior_Data(alpha_d, num_post_samples);
% [delta_mean, delta_samples] = md_post_sampling.Posterior_Discrepancy_Samples(Z);

% Obtain Optimal Solution Update via Continuation
num_continuation_steps = 3;
md_cont_update = MD_Continuation_Update(md_post_sampling, md_hessian_analysis, num_continuation_steps);
[u_cont, z_cont] = md_cont_update.Posterior_Update_Mean_PC_beta();
z_bar = z_cont(:, end);

% Sample from Posterior of Optimal Solution
% md_update = MD_Update(md_post_sampling, md_hessian_analysis);
% [z_bar, z_update_samples] = md_update.Posterior_Update_Samples();

% Sample from Design Prior of Z
% num_prior_samples = 1;
% Z_prior_samps = md_oed.Generate_Random_Design(num_prior_samples);

% Evaluate Objectives
fprintf("\n\n");
fprintf("\nObjective Value of Hi-Fi Control: \t" + opt_hifi.Jhat(z_hifi));
fprintf("\nObjective Value of Lo-Fi Control: \t" + opt_hifi.Jhat(z_lofi));
fprintf("\nObjective Value of Updated Control: \t" + opt_hifi.Jhat(z_bar));
fprintf("\n\n");
fprintf("\nError of Lo-Fi Control: \t" + oed_z_error_fn(z_lofi));
fprintf("\nError of Updated Control: \t" + oed_z_error_fn(z_bar));

% % Additional
% fprintf("\nAdditional (at Lo-Fi objective)...")
% fprintf("\nObjective Value of Hi-Fi Control: \t" + opt_lofi.Jhat(z_hifi));
% fprintf("\nObjective Value of Lo-Fi Control: \t" + opt_lofi.Jhat(z_lofi));
% fprintf("\nComparison S(z_lofi): " + abs(opt_hifi.Jhat(z_lofi)/opt_hifi.Jhat(z_hifi)-1))
% fprintf("\nComparison S_tilde(z_hifi): " + abs(opt_lofi.Jhat(z_hifi)/opt_lofi.Jhat(z_lofi)-1))

% Comparison of Lo-Fi and Hi-Fi Control Solutions
% figure;
% hold on;
% plot(x, Z_prior_samps, "Color", [0.9, 0.9, 0.7], "HandleVisibility", "off");
% plot(x, z_update_samples, "Color", [0.7, 0.9, 0.9], "HandleVisibility", "off");
% plot(x, z_lofi, "Color", 0.5 * [0.9, 0.9, 0.3], "DisplayName", "Lo-Fi Sol.");
% plot(x, z_hifi, "DisplayName", "Hi-Fi Sol.");
% plot(x, z_update_mean, "Color", 0.5 * [0.3, 0.9, 0.9], "DisplayName", "Updated Sol.");
% title("Lo-Fi & Hi-Fi Controls");
% legend("Location", "best");
% % z_prior_samples_w = z_prior_interface.Sample_with_sCovariance_W_z_Inverse(10);
% % figure;
% % plot(x, z_prior_samples)
% Comparison of Lo-Fi and Hi-Fi State Solutions
% figure;
% hold on;
% plot(x, con_hifi.State_Solve(z_lofi), "r-", "DisplayName", "Lo-Fi Sol.");
% plot(x, con_hifi.State_Solve(z_hifi), "k--", "DisplayName", "Hi-Fi Sol.");
% plot(x, con_hifi.State_Solve(z_update_mean), "b-", "DisplayName", "Updated Sol.");
% % plot(x, obj.T, "DisplayName", "Target");
% title("States for Lo-Fi & Hi-Fi Controls");
% legend("Location", "best");

% figure;
% hold on;
% plot(x, z_lofi, "r-", "DisplayName", "Lo-Fi Sol.");
% plot(x, z_hifi, "k--", "DisplayName", "Hi-Fi Sol.");
% plot(x, z_update_mean, "b-", "DisplayName", "Updated Sol.");
% % plot(x, obj.T, "DisplayName", "Target");
% title("Lo-Fi & Hi-Fi Controls");
% legend("Location", "best");
