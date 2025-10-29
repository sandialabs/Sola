%%
clear;
close all;
addpath(genpath('../../../src'));
rng(121234);

random_numbers = randn(10^6, 1);
writematrix(random_numbers, 'random_numbers.txt');

rng(121234);

m = 51;
x = linspace(0, 1, m)';

[M, S] = Assemble_Mass_and_Stiffness(m);

data_interface = MD_Data_Interface_synthetic_test_lumped_mass(m);

u_hyperparam_interface = MD_u_Hyperparameter_Interface_synthetic_test_lumped_mass(m);
u_hyperparam_interface.trace_estimator_sample_size = 100;
u_hyperparam_interface.beta_u = 0.007702351792463;
u_hyperparam_interface.alpha_d = 2.177109166165424e-07;
u_prior_interface = MD_Lumped_Mass_u_Prior_Interface(S, M, data_interface, u_hyperparam_interface);

z_hyperparam_interface = MD_z_Hyperparameter_Interface_synthetic_test_lumped_mass(m);
z_hyperparam_interface.beta_z = 0.009305846653704;
z_prior_interface = MD_Numeric_Laplacian_z_Prior_Interface(S, M, data_interface, z_hyperparam_interface, u_prior_interface);

%%
num_prior_samples = 100;
md_prior_sampling = MD_Prior_Sampling(data_interface, u_prior_interface, z_prior_interface);

%%
prior_delta_z_opt = md_prior_sampling.Prior_Discrepancy_Samples_at_z_opt(num_prior_samples);
z = zeros(m, 3);
z(:, 1) = x;
z(:, 2) = x.^2 + 1;
z(:, 3) = sin(2 * pi * x);
prior_delta = md_prior_sampling.Prior_Discrepancy_Samples(z, num_prior_samples);

%%
md_post_sampling = MD_Posterior_Sampling(data_interface, u_prior_interface, z_prior_interface);
num_post_samples = 100;
md_post_sampling.Compute_Posterior_Data(u_hyperparam_interface.alpha_d, num_post_samples);
Z_test = randn(m, 3);
Z_test(:, 1:2) = md_post_sampling.post_data.Z;
Z_test(:, 3) = 1.5 * ones(m, 1);
[post_delta_mean, post_delta_samples] = md_post_sampling.Posterior_Discrepancy_Samples(Z_test);

%%
opt_prob_interface = MD_Opt_Prob_Interface_synthetic_test_lumped_mass(m);
md_hessian_analysis = MD_Hessian_Analysis(opt_prob_interface, z_prior_interface);
num_evals = 20;
oversampling = 10;
md_hessian_analysis.Compute_Hessian_GEVP(data_interface.z_opt, num_evals, oversampling);
md_update = MD_Update(md_post_sampling, md_hessian_analysis);

[post_z_mean, post_z_samples] = md_update.Posterior_Update_Samples();

%%
post_delta_mean = reshape(cell2mat(post_delta_mean), 51, 3);
save('Sabl_Output.mat', 'prior_delta', 'prior_delta_z_opt', 'post_delta_mean', 'post_delta_samples', 'post_z_mean', 'post_z_samples');
