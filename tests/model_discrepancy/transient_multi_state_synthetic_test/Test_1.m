%%
clear;
close all;
addpath(genpath('../../../src'));
rng(121234);

suppress_figures = true;

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
for k = 1:2
    u_hyperparam_interface_cell{k} = MD_u_Hyperparameter_Interface_transient_multi_state_synthetic(n_y, n_t, k);
    u_spatial_prior_interface_cell{k} = MD_Numeric_Laplacian_u_Prior_Interface(S, M, data_interface, u_hyperparam_interface_cell{k});
    transient_prior_cov{k} = MD_Transient_Prior_Covariance_Sabl(data_interface, u_hyperparam_interface_cell{k}, 1, n_t, n_y);
    u_prior_interface_cell{k} = MD_Transient_Elliptic_u_Prior_Interface(data_interface, u_spatial_prior_interface_cell{k}, transient_prior_cov{k});
end
u_prior_interface = MD_Multi_State_u_Prior_Interface(data_interface, u_prior_interface_cell, u_hyperparam_interface_cell);

num_state_solves = 100;
z_hyperparam_interface = MD_z_Hyperparameter_Interface_transient_multi_state_synthetic(data_interface, num_state_solves, n_y, n_t);
z_prior_interface = MD_Numeric_Laplacian_z_Prior_Interface(S, M, data_interface, z_hyperparam_interface, u_prior_interface);

%%
num_prior_samples = 100;
md_prior_sampling = MD_Prior_Sampling(data_interface, u_prior_interface, z_prior_interface);
md_prior_sampling.Generate_Prior_Discrepancy_z_opt_Sample_Data(num_prior_samples);
md_prior_sampling.Generate_Prior_Discrepancy_z_pert_Sample_Data();

%%
md_post_sampling = MD_Posterior_Sampling(data_interface, u_prior_interface, z_prior_interface);
alpha_d = mean([u_hyperparam_interface_cell{1}.alpha_d; u_hyperparam_interface_cell{2}.alpha_d]);
num_post_samples = 100;
md_post_sampling.Compute_Posterior_Data(alpha_d, num_post_samples);
Z_test = randn(n_y, 3);
Z_test(:, 1:2) = md_post_sampling.post_data.Z;
Z_test(:, 3) = 1.5 * ones(n_y, 1);
[delta_mean, delta_samples] = md_post_sampling.Posterior_Discrepancy_Samples(Z_test);

%%
opt_prob_interface = MD_Opt_Prob_Interface_transient_multi_state_synthetic_test(n_y, n_t, x, M, c_low);
md_hessian_analysis = MD_Hessian_Analysis(opt_prob_interface, z_prior_interface);
md_update = MD_Update(md_post_sampling, md_hessian_analysis);

[z_update_mean, z_update_samples] = md_update.Posterior_Update_Samples();

if ~suppress_figures
    figure;
    hold on;
    for k = 1:num_post_samples
        plot(x, z_update_samples(:, k), 'color', [.9, .9, .9], 'LineWidth', 3);
    end
    plot(x, ((c_low / c_high)^(n_t / 3)) * (1 + x), 'color', 'black', 'LineWidth', 3);
    plot(x, 1 + x, 'color', 'cyan', 'LineWidth', 3);
    plot(x, z_update_mean, '--', 'color', 'red', 'LineWidth', 3);
end

%%
z_mean_ref = load('reference_solution.mat').z_update_mean;
z_samples_ref = load('reference_solution.mat').z_update_samples;
ref_diff = max(norm(z_mean_ref - z_update_mean) / norm(z_update_mean), norm(z_update_samples - z_samples_ref) / norm(z_update_samples));

if ref_diff > 1.e-7
    fprintf(2,'\nmodel_discrepancy/transient_multi_state_synthetic_test failed.\n');
else
    fprintf(1,'\nmodel_discrepancy/transient_multi_state_synthetic_test passed.\n');
end