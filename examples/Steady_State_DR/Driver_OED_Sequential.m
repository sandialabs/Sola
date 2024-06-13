% Clear Workspace and Add Interfaces to Path
clear;
close all;
% clc;
addpath(genpath('../../src'));
rng(0);

% Set Default Font Axes and Line Width
set(0, "DefaultAxesFontSize", 20);
set(0, "DefaultLineLineWidth", 3);
set(0, "DefaultLineMarkerSize", 20);

% Retrieve Model Parameters (D, Z, diff/reg/react_coeff, m, u_lofi, z_hifi/lofi; remove Z and D though)
load Optimization_Results.mat;
clear Z D;

% Set Hi-Fi and Lo-Fi Objectives and Constraints
obj = Diff_React_Objective(m, reg_coeff);
con_lofi = Diff_React_Constraint(m, diff_coeff, react_coeff);
opt_lofi = Reduced_Space_Optimization(obj, con_lofi);
con_hifi = Diff_React_HiFi_Constraint(con_lofi);
opt_hifi = Reduced_Space_Optimization(obj, con_hifi);
x = con_lofi.x;

best_objective = opt_hifi.Jhat(z_hifi);
% Show initial objective
fprintf("\nStep 0:\n-------------");
new_objective = opt_hifi.Jhat(z_lofi);
fprintf('Objective of z_lofi: \t%.2f\n', new_objective);
fprintf('Objective of z_hifi: \t%.2f\n\n', best_objective);
old_objective = new_objective;

% Set Data Interface (no data there yet, except for z_lofi/u_lofi)
data_interface = MD_Data_Interface_Diff_React(u_lofi, z_lofi);

% Generate Priors for u and z
alpha_u = 2^2;
alpha_z = 1.e-10;
u_prior_interface = MD_Elliptic_u_Prior_Interface_Diff_React(alpha_u, opt_lofi);
z_prior_interface = MD_Elliptic_z_Prior_Interface_Diff_React(alpha_z, opt_lofi);

% Error with z_hifi
oed_z_error_fn = @(z) sqrt((z - z_hifi)' * z_prior_interface.Apply_M_z(z - z_hifi)) / sqrt(z_hifi' * z_prior_interface.Apply_M_z(z_hifi));

% %% Hessian analysis
opt_prob_interface = MD_Opt_Prob_Interface_Sabl(opt_lofi, data_interface);
md_hessian_analysis = MD_Hessian_Analysis(opt_prob_interface, z_prior_interface);
num_evals = 4;
oversampling = 20;
md_hessian_analysis.Compute_Hessian_GEVP(data_interface.z_opt, num_evals, oversampling);

% Perform Offline OED Computations - This is used to generate many random designs (to avoid OED in next steps)
alpha_zd = 1.e-2;
beta_zd = 1.e-2;
oed_interface = MD_OED_Interface_Diff_React(data_interface, con_lofi, alpha_zd, beta_zd);
md_oed = MD_OED(opt_prob_interface, data_interface, u_prior_interface, z_prior_interface, md_hessian_analysis, oed_interface);
md_oed.Offline_Computation();
Z_random_samples = md_oed.Generate_Random_Design(100);

% Base Case
fprintf("\nStep 1:\n-------------");
Z = Z_random_samples(:, 1);
D = Evaluate_Discrepancy(con_hifi, con_lofi, Z);
data_interface.Set_Z_and_D(Z, D);
alpha_d = 1.e-4;

% Obtain Optimal Solution Update
md_post_sampling = MD_Posterior_Sampling(data_interface, u_prior_interface, z_prior_interface);
md_post_sampling.Compute_Posterior_Data(alpha_d, 1);
md_update = MD_Update(md_post_sampling, md_hessian_analysis);
z_update_mean = md_update.Posterior_Update_Samples();

% Display Stats
new_objective = opt_hifi.Jhat(z_update_mean);
fprintf('Objective of z_lofi: \t%.2f\n', new_objective);
fprintf('Percent Improvement: \t%.2f%%\n\n', 100 * (old_objective - new_objective) / (old_objective - best_objective));
old_objective = new_objective;

%% STEP 2 & Beyond
N = 6;
for p = 2:N
    % Update Data Interface (BEWARE: the change affects the dependencies directly!)
    fprintf('\nStep %d:\n-------------\n', p);
    updated_data_interface = MD_Data_Interface_Diff_React(u_lofi, z_update_mean);

    % Obtain Discrepancies
    Z = Z_random_samples(:, 1:p);
    D = Evaluate_Discrepancy(con_hifi, con_lofi, Z);
    updated_data_interface.Set_Z_and_D(Z, D);

    % % Obtain Optimal Solution Update
    md_post_sampling = MD_Posterior_Sampling(updated_data_interface, u_prior_interface, z_prior_interface);
    md_post_sampling.Compute_Posterior_Data(alpha_d, 1);
    updated_data_interface.z_opt = z_lofi;
    updated_data_interface.u_opt = u_lofi;
    md_update = MD_Update(md_post_sampling, md_hessian_analysis);
    z_update_mean = md_update.Posterior_Update_Samples();

    % Display Stats
    new_objective = opt_hifi.Jhat(z_update_mean);
    fprintf('Objective of z_bar: \t%.2f\n', new_objective);
    fprintf('Rel. Err of z_bar: \t%.2f%%\n\n', 100 * oed_z_error_fn(z_update_mean));
    fprintf('Percent Improvement: \t%.2f%%\n\n', 100 * (old_objective - new_objective) / (old_objective - best_objective));
    old_objective = new_objective;
end
