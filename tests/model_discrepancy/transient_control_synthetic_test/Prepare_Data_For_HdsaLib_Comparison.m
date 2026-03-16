%%
clear;
close all;
rng(121234);

random_numbers = randn(10^6, 1);
writematrix(random_numbers, 'random_numbers.txt');

rng(121234);

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

%%
md_post_sampling = MD_Posterior_Sampling(data_interface, u_prior_interface, z_prior_interface);
alpha_d = u_hyperparam_interface.alpha_d;
num_post_samples = 100;
md_post_sampling.Compute_Posterior_Data(alpha_d, num_post_samples);

%%
md_hessian_analysis = MD_Hessian_Analysis(opt_prob_interface, z_prior_interface);
num_evals = 30;
oversampling = 9;
md_hessian_analysis.Compute_Hessian_GEVP(data_interface.z_opt, num_evals, oversampling);
md_update = MD_Update(md_post_sampling, md_hessian_analysis);

[post_z_mean, post_z_samples] = md_update.Posterior_Update_Samples();

%%
save('Sabl_Output.mat', 'post_z_mean', 'post_z_samples');
