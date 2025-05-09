%%
clear;
close all;
addpath(genpath('../../../src'));
rng(121234);

random_numbers = randn(10^6, 1);
writematrix(random_numbers, 'random_numbers.txt');

rng(121234);

n_y = 50;
n_t = 10;
[M, S, x] = Assemble_Mass_and_Stiffness(n_y);
c_low = 0.95;
c_high = 0.93;

data_interface = MD_Data_Interface_transient_multi_state_synthetic(n_y, n_t, c_low, c_high);

u_hyperparam_interface_cell = cell(2, 1);
u_spatial_prior_interface_cell = cell(2, 1);
transient_prior_cov = cell(2, 1);
u_prior_interface_cell = cell(2, 1);

u_hyperparam_interface_cell{1} = MD_u_Hyperparameter_Interface_transient_multi_state_synthetic(n_y, n_t, 1);
%u_hyperparam_interface_cell{1}.alpha_u = 0.001688110759857;
u_hyperparam_interface_cell{1}.beta_u = 0.009166435191031;
u_hyperparam_interface_cell{1}.beta_t = 0.027499305573092;
u_spatial_prior_interface_cell{1} = MD_Numeric_Laplacian_u_Prior_Interface(S, M, data_interface, u_hyperparam_interface_cell{1});
transient_prior_cov{1} = MD_Transient_Prior_Covariance_Sabl(data_interface, u_hyperparam_interface_cell{1}, 1, n_t, n_y);
u_prior_interface_cell{1} = MD_Transient_Elliptic_u_Prior_Interface(data_interface, u_spatial_prior_interface_cell{1}, transient_prior_cov{1});

u_hyperparam_interface_cell{2} = MD_u_Hyperparameter_Interface_transient_multi_state_synthetic(n_y, n_t, 2);
%u_hyperparam_interface_cell{2}.alpha_u = 0.006235002943316;
u_hyperparam_interface_cell{2}.beta_u = 0.009166435191031;
u_hyperparam_interface_cell{2}.beta_t = 0.027499305573092;
u_spatial_prior_interface_cell{2} = MD_Numeric_Laplacian_u_Prior_Interface(S, M, data_interface, u_hyperparam_interface_cell{2});
transient_prior_cov{2} = MD_Transient_Prior_Covariance_Sabl(data_interface, u_hyperparam_interface_cell{2}, 1, n_t, n_y);
u_prior_interface_cell{2} = MD_Transient_Elliptic_u_Prior_Interface(data_interface, u_spatial_prior_interface_cell{2}, transient_prior_cov{2});

u_prior_interface = MD_Multi_State_u_Prior_Interface(data_interface, u_prior_interface_cell, u_hyperparam_interface_cell);

num_state_solves = 100;
z_hyperparam_interface = MD_z_Hyperparameter_Interface_transient_multi_state_synthetic(data_interface, num_state_solves, n_y, n_t);
%z_hyperparam_interface.alpha_z = 1.076021648025798e+03;
z_hyperparam_interface.beta_z = 0.009305846653704;
z_prior_interface = MD_Numeric_Laplacian_z_Prior_Interface(S, M, data_interface, z_hyperparam_interface, u_prior_interface);

%%
num_prior_samples = 100;
md_prior_sampling = MD_Prior_Sampling(data_interface, u_prior_interface, z_prior_interface);

prior_delta_z_opt = md_prior_sampling.Prior_Discrepancy_Samples_at_z_opt(num_prior_samples);

z = zeros(n_y, 3);
z(:, 1) = x;
z(:, 2) = x.^2 + 1;
z(:, 3) = sin(2 * pi * x);
prior_delta = md_prior_sampling.Prior_Discrepancy_Samples(z, num_prior_samples);

%%
md_post_sampling = MD_Posterior_Sampling(data_interface, u_prior_interface, z_prior_interface);
alpha_d = mean([u_hyperparam_interface_cell{1}.alpha_d; u_hyperparam_interface_cell{2}.alpha_d]);
num_post_samples = 100;
md_post_sampling.Compute_Posterior_Data(alpha_d, num_post_samples);
Z_test = randn(n_y, 3);
Z_test(:, 1:2) = md_post_sampling.post_data.Z;
Z_test(:, 3) = 1.5 * ones(n_y, 1);
[post_delta_mean, post_delta_samples] = md_post_sampling.Posterior_Discrepancy_Samples(Z_test);

%%
opt_prob_interface = MD_Opt_Prob_Interface_transient_multi_state_synthetic_test(n_y, n_t, x, M, c_low);
md_hessian_analysis = MD_Hessian_Analysis(opt_prob_interface, z_prior_interface);
md_update = MD_Update(md_post_sampling, md_hessian_analysis);

[post_z_mean, post_z_samples] = md_update.Posterior_Update_Samples();

%%
post_delta_mean = reshape(cell2mat(post_delta_mean), 1000, 3);
save('Sabl_Output.mat', 'prior_delta', 'prior_delta_z_opt', 'post_delta_mean', 'post_delta_samples', 'post_z_mean', 'post_z_samples');
