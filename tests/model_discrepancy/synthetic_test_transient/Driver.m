%%
clear;
close all;
addpath(genpath('../../../src'));
rng(121234);

suppress_figures = true;

n_y = 51;
n_t = 10;
c_low = 0.95;
c_high = 0.93;
T = 1;
x = linspace(0, 1, n_y)';

[M, S] = Assemble_Mass_and_Stiffness(n_y);

data_interface = MD_Data_Interface_synthetic_test_transient(n_y, n_t, T, c_low, c_high);

u_hyperparam_interface = MD_u_Hyperparameter_Interface_synthetic_test_transient(n_y, n_t, T);
spatial_u_prior_interface = MD_Numeric_Laplacian_u_Prior_Interface(S, M, data_interface, u_hyperparam_interface);
transient_prior_cov = MD_Transient_Prior_Covariance_Sabl(data_interface, u_hyperparam_interface, T, n_t, n_y);
u_prior_interface = MD_Transient_Elliptic_u_Prior_Interface(data_interface, spatial_u_prior_interface, transient_prior_cov);

z_hyperparam_interface = MD_z_Hyperparameter_Interface_synthetic_test_transient(n_y);
z_prior_interface = MD_Numeric_Laplacian_z_Prior_Interface(S, M, data_interface, z_hyperparam_interface, u_prior_interface);

opt_prob_interface = MD_Opt_Prob_Interface_synthetic_test_transient(n_y, n_t, T, c_low);

%%
num_prior_samples = 100;
md_prior_sampling = MD_Prior_Sampling(data_interface, u_prior_interface, z_prior_interface);

%%
delta_prior_samples_zopt = md_prior_sampling.Prior_Discrepancy_Samples_at_z_opt(num_prior_samples);

z = zeros(n_y, 3);
z(:, 1) = x;
z(:, 2) = x.^2 + 1;
z(:, 3) = sin(2 * pi * x);
delta_prior_samples = md_prior_sampling.Prior_Discrepancy_Samples(z, num_prior_samples);
%%
md_post_sampling = MD_Posterior_Sampling(data_interface, u_prior_interface, z_prior_interface);
num_post_samples = 100;
md_post_sampling.Compute_Posterior_Data(u_hyperparam_interface.alpha_d, num_post_samples);
Z_test = randn(n_y, 3);
Z_test(:, 1:2) = md_post_sampling.post_data.Z;
Z_test(:, 3) = 1.5 * ones(n_y, 1);
[delta_mean, delta_samples] = md_post_sampling.Posterior_Discrepancy_Samples(Z_test);

if ~suppress_figures
    k = 2;
    ts = 8;

    d_mean = reshape(delta_mean{k}, n_y, n_t);
    d_samples = reshape(delta_samples{k}, n_y, n_t, num_post_samples);
    d_data = reshape(data_interface.D(:, k), n_y, n_t);

    figure;
    hold on;
    plot(x, reshape(d_samples(:, ts, :), n_y, num_post_samples), 'color', [.9, .9, .9], 'LineWidth', 3);
    plot(x, d_data(:, ts), 'color', 'black', 'LineWidth', 3);
    plot(x, d_mean(:, ts), '--', 'Color', 'red', 'LineWidth', 3);
end

%%
md_hessian_analysis = MD_Hessian_Analysis(opt_prob_interface, z_prior_interface);
num_evals = 20;
oversampling = 10;
md_hessian_analysis.Compute_Hessian_GEVP(data_interface.z_opt, num_evals, oversampling);
md_update = MD_Update(md_post_sampling, md_hessian_analysis);

[z_update_mean, z_update_samples] = md_update.Posterior_Update_Samples();

if ~suppress_figures
    figure;
    hold on;
    plot(x, (c_low / c_high)^((n_t - 1) / 3) * (1 + x), 'color', 'black', 'LineWidth', 3);
    plot(x, 1 + x, 'color', 'cyan', 'LineWidth', 3);
    plot(x, z_update_mean, '--', 'color', 'red', 'LineWidth', 3);
    for k = 1:num_post_samples
        plot(x, z_update_samples(:, k), 'color', [.9, .9, .9], 'LineWidth', 3);
    end
    plot(x, (c_low / c_high)^((n_t - 1) / 3) * (1 + x), 'color', 'black', 'LineWidth', 3);
    plot(x, 1 + x, 'color', 'cyan', 'LineWidth', 3);
    plot(x, z_update_mean, '--', 'color', 'red', 'LineWidth', 3);
end

%%
z_mean_ref = load('reference_solution.mat').z_update_mean;
z_samples_ref = load('reference_solution.mat').z_update_samples;
ref_diff = max(norm(z_mean_ref - z_update_mean) / norm(z_update_mean), norm(z_update_samples - z_samples_ref) / norm(z_update_samples));

if ref_diff > 1.e-9
    fprintf(2,'\nModel discrepancy synthetic_test_transient failed.\n');
else
    fprintf(1,'\nModel discrepancy synthetic_test_transient passed.\n');
end