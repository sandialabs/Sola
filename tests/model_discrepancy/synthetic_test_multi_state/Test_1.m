%%
clear;
close all;
rng(121234);

suppress_figures = true;

m = 51;
[M, S, x] = Assemble_Mass_and_Stiffness(m);
state_map_array = cell(2, 1);
state_map_array{1} = (1:m)';
state_map_array{2} = ((m + 1):(2 * m))';

data_interface = MD_Data_Interface_multi_state_synthetic_test(m);

u_hyperparam_interface_cell = cell(2, 1);
u_prior_interface_cell = cell(2, 1);

u_hyperparam_interface_cell{1} = MD_u_Hyperparameter_Interface_multi_state_synthetic_test(1, m);
u_prior_interface_cell{1} = MD_Numeric_Laplacian_u_Prior_Interface(S, M, data_interface, u_hyperparam_interface_cell{1});

u_hyperparam_interface_cell{2} = MD_u_Hyperparameter_Interface_multi_state_synthetic_test(2, m);
u_prior_interface_cell{2} = MD_Numeric_Laplacian_u_Prior_Interface(S, M, data_interface, u_hyperparam_interface_cell{2});

u_prior_interface = MD_Multi_State_u_Prior_Interface(data_interface, u_prior_interface_cell, u_hyperparam_interface_cell);

z_hyperparam_interface = MD_z_Hyperparameter_Interface_multi_state_synthetic_test(m);
z_prior_interface = MD_Numeric_Laplacian_z_Prior_Interface(S, M, data_interface, z_hyperparam_interface, u_prior_interface);

%%
num_prior_samples = 100;
md_prior_sampling = MD_Prior_Sampling(data_interface, u_prior_interface, z_prior_interface);
md_prior_sampling.Generate_Prior_Discrepancy_Sample_Data(num_prior_samples);

%%
md_post_sampling = MD_Posterior_Sampling(data_interface, u_prior_interface, z_prior_interface);
alpha_d = mean([u_hyperparam_interface_cell{1}.alpha_d(); u_hyperparam_interface_cell{2}.alpha_d()]);
num_post_samples = 100;
md_post_sampling.Compute_Posterior_Data(alpha_d, num_post_samples);
Z_test = randn(m, 3);
Z_test(:, 1:2) = md_post_sampling.post_data.Z;
Z_test(:, 3) = 1.5 * ones(m, 1);
[delta_mean, delta_samples] = md_post_sampling.Posterior_Discrepancy_Samples(Z_test);

if ~suppress_figures

    figure;
    subplot(2, 1, 1);
    hold on;
    for k = 1:num_post_samples
        plot(x, delta_samples{1}(1:m, k), 'color', [.9, .9, .9], 'LineWidth', 3);
    end
    plot(x, md_post_sampling.post_data.D(1:m, 1), 'color', 'black', 'LineWidth', 3);
    plot(x, delta_mean{1}(1:m), '--', 'color', 'red', 'LineWidth', 3);
    subplot(2, 1, 2);
    hold on;
    for k = 1:num_post_samples
        plot(x, delta_samples{1}((m + 1):end, k), 'color', [.9, .9, .9], 'LineWidth', 3);
    end
    plot(x, md_post_sampling.post_data.D((m + 1):end, 1), 'color', 'black', 'LineWidth', 3);
    plot(x, delta_mean{1}((m + 1):end), '--', 'color', 'red', 'LineWidth', 3);

    figure;
    subplot(2, 1, 1);
    hold on;
    for k = 1:num_post_samples
        plot(x, delta_samples{2}(1:m, k), 'color', [.9, .9, .9], 'LineWidth', 3);
    end
    plot(x, md_post_sampling.post_data.D(1:m, 2), 'color', 'black', 'LineWidth', 3);
    plot(x, delta_mean{2}(1:m), '--', 'color', 'red', 'LineWidth', 3);
    subplot(2, 1, 2);
    hold on;
    for k = 1:num_post_samples
        plot(x, delta_samples{2}((m + 1):end, k), 'color', [.9, .9, .9], 'LineWidth', 3);
    end
    plot(x, md_post_sampling.post_data.D((m + 1):end, 2), 'color', 'black', 'LineWidth', 3);
    plot(x, delta_mean{2}((m + 1):end), '--', 'color', 'red', 'LineWidth', 3);

    figure;
    subplot(2, 1, 1);
    hold on;
    for k = 1:num_post_samples
        plot(x, delta_samples{3}(1:m, k), 'color', [.9, .9, .9], 'LineWidth', 3);
    end
    plot(x, delta_mean{3}(1:m), '--', 'color', 'red', 'LineWidth', 3);
    subplot(2, 1, 2);
    hold on;
    for k = 1:num_post_samples
        plot(x, delta_samples{3}((m + 1):end, k), 'color', [.9, .9, .9], 'LineWidth', 3);
    end
    plot(x, delta_mean{3}((m + 1):end), '--', 'color', 'red', 'LineWidth', 3);

end
%%
opt_prob_interface = MD_Opt_Prob_Interface_multi_state_synthetic_test(m);
md_hessian_analysis = MD_Hessian_Analysis(opt_prob_interface, z_prior_interface);
md_update = MD_Update(md_post_sampling, md_hessian_analysis);

[z_update_mean, z_update_samples] = md_update.Posterior_Update_Samples();

if ~suppress_figures
    figure;
    hold on;
    for k = 1:num_post_samples
        plot(x, z_update_samples(:, k), 'color', [.9, .9, .9], 'LineWidth', 3);
    end
    plot(x, (1 + x) / (1.1^(2 / 3)), 'color', 'black', 'LineWidth', 3);
    plot(x, 1 + x, 'color', 'cyan', 'LineWidth', 3);
    plot(x, z_update_mean, '--', 'color', 'red', 'LineWidth', 3);
end

%%
z_mean_ref = load('reference_solution.mat').z_update_mean;
z_samples_ref = load('reference_solution.mat').z_update_samples;
ref_diff = max(norm(z_mean_ref - z_update_mean) / norm(z_update_mean), norm(z_update_samples - z_samples_ref) / norm(z_update_samples));

if ref_diff > 1.e-9
    fprintf(2, '\nmodel_discrepancy/synthetic_test_multi_state failed.\n');
else
    fprintf(1, '\nmodel_discrepancy/synthetic_test_multi_state passed.\n');
end
