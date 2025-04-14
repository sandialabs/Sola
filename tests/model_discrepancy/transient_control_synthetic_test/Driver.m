%%
clear;
close all;
clc;
addpath(genpath('../../../src'));
rng(121234);

suppress_figures = true;

n_y = 50;
x = linspace(0, 1, n_y)';
n_t = 20;
T = 1;
t = linspace(0, T, n_t)';

data_interface = MD_Data_Interface_transient_control_synthetic_test(n_y, n_t);
opt_prob_interface = MD_Opt_Prob_Interface_transient_control_synthetic_test(n_y, n_t);

[M, S] = Assemble_Mass_and_Stiffness(n_y);
u_hyperparam_interface = MD_u_Hyperparameter_Interface_transient_control_synthetic_test(n_y, n_t);
transient_prior_cov = MD_Transient_Prior_Covariance_Sabl(data_interface, u_hyperparam_interface, T, n_t, n_y);
spatial_u_prior_interface = MD_Numeric_Laplacian_u_Prior_Interface(S, M, data_interface, u_hyperparam_interface);
u_prior_interface = MD_Transient_Elliptic_u_Prior_Interface(data_interface, spatial_u_prior_interface, transient_prior_cov);

[Mt, St] = Assemble_Mass_and_Stiffness(n_t);
num_controls = 2;
num_state_solves = 100;
z_hyperparam_interface = MD_z_Hyperparameter_Interface_transient_control_synthetic_test(num_state_solves, t, opt_prob_interface);
z_prior_interface = MD_Transient_Vector_z_Prior_Interface(St, Mt, num_controls, data_interface, z_hyperparam_interface, u_prior_interface);

%%
num_prior_samples = 100;
md_prior_sampling = MD_Prior_Sampling(data_interface, u_prior_interface, z_prior_interface);
md_prior_sampling.Generate_Prior_Discrepancy_Sample_Data(num_prior_samples);

%md_prior_vis = MD_Prior_Visualization(md_prior_sampling);
%md_prior_vis.Visualization_for_Prior_Time_Evolution(1, true);
%md_prior_vis.Visualization_for_Prior_Discrepancy_at_z_opt(1);
%md_prior_vis.Visualization_for_Prior_Discrepancy_at_z_pert(1);

%%
md_post_sampling = MD_Posterior_Sampling(data_interface, u_prior_interface, z_prior_interface);
alpha_d = u_hyperparam_interface.alpha_d;
num_post_samples = 100;
md_post_sampling.Compute_Posterior_Data(alpha_d, num_post_samples);

%%
md_hessian_analysis = MD_Hessian_Analysis(opt_prob_interface, z_prior_interface);
num_evals = 30;
oversampling = 9;
md_hessian_analysis.Compute_Hessian_GEVP(data_interface.z_opt,num_evals,oversampling);
md_update = MD_Update(md_post_sampling, md_hessian_analysis);

[z_update_mean, z_update_samples] = md_update.Posterior_Update_Samples();

%%
% save('reference_solution.mat','z_hyperparam_interface','md_hessian_analysis')
ref = load('reference_solution.mat');
ref_diff = norm(z_hyperparam_interface.alpha_z-ref.z_hyperparam_interface.alpha_z);
ref_diff = min(ref_diff,norm(z_hyperparam_interface.beta_t-ref.z_hyperparam_interface.beta_t));
ref_diff = min(ref_diff,norm(md_hessian_analysis.evals-ref.md_hessian_analysis.evals));

if ref_diff > 1.e-9
    disp('model_discrepancy_sythetic_test difference:');
    disp(ref_diff);
end
