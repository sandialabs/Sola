%%
clear;
close all;
addpath(genpath('../../../src'));
rng(121234);

suppress_figures = true;

m = 51;
x = linspace(0, 1, m)';

[M, S] = Assemble_Mass_and_Stiffness(m);

data_interface = MD_Data_Interface_synthetic_test_with_hyperparam(m);

u_hyperparam_interface = MD_u_Hyperparameter_Interface_synthetic_test_with_hyperparam(m);
u_prior_interface = MD_Numeric_Laplacian_u_Prior_Interface(S, M, data_interface, u_hyperparam_interface);

z_hyperparam_interface = MD_z_Hyperparameter_Interface_synthetic_test_with_hyperparam(m);
z_prior_interface = MD_Numeric_Laplacian_z_Prior_Interface(S, M, data_interface, z_hyperparam_interface, u_prior_interface);

%%
num_prior_samples = 100;
md_prior_sampling = MD_Prior_Sampling(data_interface, u_prior_interface, z_prior_interface);

%%
delta_prior_samples_zopt = md_prior_sampling.Prior_Discrepancy_Samples_at_z_opt(num_prior_samples);

if ~suppress_figures
    figure;
    plot(x, delta_prior_samples_zopt, 'LineWidth', 3, 'color', [.9, .9, .9]);
end

z = zeros(m, 3);
z(:, 1) = x;
z(:, 2) = x.^2 + 1;
z(:, 3) = sin(2 * pi * x);
delta_prior_samples = md_prior_sampling.Prior_Discrepancy_Samples(z, num_prior_samples);
if ~suppress_figures
    for k = 1:10
        figure;
        hold on;
        plot(x, delta_prior_samples{k}, 'LineWidth', 3);
    end
end

%%
md_post_sampling = MD_Posterior_Sampling(data_interface, u_prior_interface, z_prior_interface);
num_post_samples = 100;
md_post_sampling.Compute_Posterior_Data(u_hyperparam_interface.alpha_d, num_post_samples);
Z_test = randn(m, 3);
Z_test(:, 1:2) = md_post_sampling.post_data.Z;
Z_test(:, 3) = 1.5 * ones(m, 1);
[delta_mean, delta_samples] = md_post_sampling.Posterior_Discrepancy_Samples(Z_test);

if ~suppress_figures
    figure;
    hold on;
    plot(x, md_post_sampling.post_data.D(:, 1), 'color', 'black', 'LineWidth', 3);
    plot(x, delta_mean{1}, '--', 'color', 'red', 'LineWidth', 3);
    for k = 1:num_post_samples
        plot(x, delta_samples{1}(:, k), 'color', [.9, .9, .9], 'LineWidth', 3);
    end
    plot(x, md_post_sampling.post_data.D(:, 1), 'color', 'black', 'LineWidth', 3);
    plot(x, delta_mean{1}, '--', 'color', 'red', 'LineWidth', 3);

    figure;
    hold on;
    plot(x, md_post_sampling.post_data.D(:, 2), 'color', 'black', 'LineWidth', 3);
    plot(x, delta_mean{2}, '--', 'color', 'red', 'LineWidth', 3);
    for k = 1:num_post_samples
        plot(x, delta_samples{2}(:, k), 'color', [.9, .9, .9], 'LineWidth', 3);
    end
    plot(x, md_post_sampling.post_data.D(:, 2), 'color', 'black', 'LineWidth', 3);
    plot(x, delta_mean{2}, '--', 'color', 'red', 'LineWidth', 3);

    figure;
    hold on;
    plot(x, delta_mean{3}, '--', 'color', 'red', 'LineWidth', 3);
    for k = 1:num_post_samples
        plot(x, delta_samples{3}(:, k), 'color', [.9, .9, .9], 'LineWidth', 3);
    end
    plot(x, delta_mean{3}, '--', 'color', 'red', 'LineWidth', 3);

end

%%
opt_prob_interface = MD_Opt_Prob_Interface_synthetic_test_with_hyperparam(m);
md_hessian_analysis = MD_Hessian_Analysis(opt_prob_interface, z_prior_interface);
num_evals = 20;
oversampling = 10;
md_hessian_analysis.Compute_Hessian_GEVP(data_interface.z_opt, num_evals, oversampling);
md_update = MD_Update(md_post_sampling, md_hessian_analysis);

[z_update_mean, z_update_samples] = md_update.Posterior_Update_Samples();

if ~suppress_figures
    figure;
    hold on;
    plot(x, (1 + x) / (1.2^(1 / 3)), 'color', 'black', 'LineWidth', 3);
    plot(x, 1 + x, 'color', 'cyan', 'LineWidth', 3);
    plot(x, z_update_mean, '--', 'color', 'red', 'LineWidth', 3);
    for k = 1:num_post_samples
        plot(x, z_update_samples(:, k), 'color', [.9, .9, .9], 'LineWidth', 3);
    end
    plot(x, (1 + x) / (1.2^(1 / 3)), 'color', 'black', 'LineWidth', 3);
    plot(x, 1 + x, 'color', 'cyan', 'LineWidth', 3);
    plot(x, z_update_mean, '--', 'color', 'red', 'LineWidth', 3);
end

%%
z_mean_ref = load('reference_solution.mat').z_update_mean;
z_samples_ref = load('reference_solution.mat').z_update_samples;
ref_diff = max(norm(z_mean_ref - z_update_mean) / norm(z_update_mean), norm(z_update_samples - z_samples_ref) / norm(z_update_samples));
if ref_diff > 1.e-9
    disp('model_discrepancy_sythetic_test_with_hyperparam_auto_1D:');
    disp(ref_diff);
end
