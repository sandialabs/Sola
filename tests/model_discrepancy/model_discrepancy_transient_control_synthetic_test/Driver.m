%%
clear;
close all;
clc
addpath(genpath('../../../src'));
rng(121234);

suppress_figures = true;

n_y = 50;
x = linspace(0, 1, n_y)';
n_t = 20;
T = 1;
t = linspace(0,T,n_t)';

data_interface = MD_Data_Interface_transient_control_synthetic_test(n_y,n_t);
opt_prob_interface = MD_Opt_Prob_Interface_transient_control_synthetic_test(n_y,n_t);

[M, S] = Assemble_Mass_and_Stiffness(n_y);
u_hyperparams = MD_u_Hyperparameters_transient_control_synthetic_test(data_interface, n_y);
transient_prior_cov = MD_Transient_Prior_Covariance_Sabl(u_hyperparams, T, n_t, n_y);
spatial_u_prior_interface = MD_Numeric_Laplacian_u_Prior_Interface(S, M, u_hyperparams);
u_prior_interface = MD_Transient_Elliptic_u_Prior_Interface(spatial_u_prior_interface, transient_prior_cov);

[Mt, St] = Assemble_Mass_and_Stiffness(n_t);
num_controls = 2;
num_state_solves = 100;
z_hyperparams = MD_z_Hyperparameters_transient_control_synthetic_test(data_interface, u_prior_interface, num_state_solves, t, opt_prob_interface);
z_prior_interface = MD_Transient_Vector_z_Prior_Interface(St, Mt, num_controls, z_hyperparams);

%%
num_prior_samples = 100;
md_prior_sampling = MD_Prior_Sampling(data_interface, u_prior_interface, z_prior_interface);
md_prior_sampling.Generate_Prior_Discrepancy_Sample_Data(num_prior_samples);

%md_prior_sampling.Visualization_for_Prior_Time_Evolution(1);
%md_prior_sampling.Visualization_for_Prior_Discrepancy_at_z_opt(1);
md_prior_sampling.Visualization_for_Prior_Discrepancy_at_z_pert(1);

%%
md_post_sampling = MD_Posterior_Sampling(data_interface, u_prior_interface, z_prior_interface);
u_hyperparams.Determine_alpha_d();
alpha_d = u_hyperparams.alpha_d;
num_post_samples = 100;
md_post_sampling.Compute_Posterior_Data(alpha_d, num_post_samples);

%%
md_hessian_analysis = MD_Hessian_Analysis(opt_prob_interface, z_prior_interface);
md_update = MD_Update(md_post_sampling, md_hessian_analysis);

[z_update_mean, z_update_samples] = md_update.Posterior_Update_Samples();

% if ~suppress_figures
%     figure;
%     hold on;
%     plot(x, (1 + x) / (1.2^(1 / 3)), 'color', 'black', 'LineWidth', 3);
%     plot(x, 1 + x, 'color', 'cyan', 'LineWidth', 3);
%     plot(x, z_update_mean, '--', 'color', 'red', 'LineWidth', 3);
%     for k = 1:num_post_samples
%         plot(x, z_update_samples(:, k), 'color', [.9, .9, .9], 'LineWidth', 3);
%     end
%     plot(x, (1 + x) / (1.2^(1 / 3)), 'color', 'black', 'LineWidth', 3);
%     plot(x, 1 + x, 'color', 'cyan', 'LineWidth', 3);
%     plot(x, z_update_mean, '--', 'color', 'red', 'LineWidth', 3);
% end
% 
% %%
% z_mean_ref = load('reference_solution.mat').z_update_mean;
% z_samples_ref = load('reference_solution.mat').z_update_samples;
% ref_diff = max(norm(z_mean_ref - z_update_mean) / norm(z_update_mean), norm(z_update_samples - z_samples_ref) / norm(z_update_samples));
% if ref_diff > 1.e-9
%     disp('model_discrepancy_sythetic_test difference:');
%     disp(ref_diff);
% end
