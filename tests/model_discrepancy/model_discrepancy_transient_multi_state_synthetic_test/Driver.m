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
c_high = 0.92;

data_interface = MD_Data_Interface_transient_multi_state_synthetic_test(n_y,n_t,c_low,c_high);
data_interface.Load_Data();

ui_hyperparams = cell(2,1);
ui_spatial_prior_interface = cell(2,1);
transient_prior_cov = cell(2,1);
ui_prior_interface = cell(2,1);
for k = 1:2
    ui_hyperparams{k} = MD_u_Hyperparameters_transient_multi_state_synthetic_test(data_interface, n_y, n_t, k);
    ui_spatial_prior_interface{k} = MD_Numeric_Laplacian_u_Prior_Interface(S,M,ui_hyperparams{k});
    transient_prior_cov{k} = MD_Transient_Prior_Covariance_Sabl(ui_hyperparams{k}, 1, n_t, n_y);
    ui_prior_interface{k} = MD_Transient_Elliptic_u_Prior_Interface(ui_spatial_prior_interface{k}, transient_prior_cov{k});
end
u_prior_interface = MD_Multi_State_u_Prior_Interface(ui_prior_interface);

num_state_solves = 100;
z_hyperparams = MD_z_Hyperparameters_transient_multi_state_synthetic_test(data_interface, u_prior_interface, num_state_solves, n_y, n_t);
z_prior_interface = MD_Numeric_Laplacian_z_Prior_Interface(S, M, z_hyperparams);

%%
num_prior_samples = 100;
md_prior_sampling = MD_Prior_Sampling(data_interface, u_prior_interface, z_prior_interface);

delta_samples = md_prior_sampling.Prior_Discrepancy_Samples_at_z_opt(num_prior_samples);

if ~suppress_figures
    figure;
    plot(x, delta_samples, 'LineWidth', 3, 'color', [.9, .9, .9]);
end

z = zeros(n_y, 3);
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
ui_hyperparams{1}.Determine_alpha_d();
ui_hyperparams{2}.Determine_alpha_d();
alpha_d = mean([ui_hyperparams{1}.alpha_d;ui_hyperparams{2}.alpha_d]);
num_post_samples = 100;
md_post_sampling.Compute_Posterior_Data(alpha_d, num_post_samples);
Z_test = randn(n_y, 3);
Z_test(:, 1:2) = md_post_sampling.post_data.Z;
Z_test(:, 3) = 1.5 * ones(n_y, 1);
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
opt_prob_interface = MD_Opt_Prob_Interface_transient_multi_state_synthetic_test(n_y,n_t,x,M,c_low);
md_hessian_analysis = MD_Hessian_Analysis(opt_prob_interface, z_prior_interface);
md_update = MD_Update(md_post_sampling, md_hessian_analysis);

[z_update_mean, z_update_samples] = md_update.Posterior_Update_Samples();

if true %~suppress_figures
    figure;
    hold on;
    for k = 1:num_post_samples
        plot(x, z_update_samples(:, k), 'color', [.9, .9, .9], 'LineWidth', 3);
    end
    plot(x, ( (c_low/c_high)^(n_t/3) ) * (1 + x), 'color', 'black', 'LineWidth', 3);
    plot(x, 1 + x, 'color', 'cyan', 'LineWidth', 3);
    plot(x, z_update_mean, '--', 'color', 'red', 'LineWidth', 3);
end

% %%
% z_mean_ref = load('reference_solution.mat').z_update_mean;
% z_samples_ref = load('reference_solution.mat').z_update_samples;
% ref_diff = max(norm(z_mean_ref - z_update_mean) / norm(z_update_mean), norm(z_update_samples - z_samples_ref) / norm(z_update_samples));
% if ref_diff > 1.e-9
%     disp('model_discrepancy_sythetic_test difference:');
%     disp(ref_diff);
% end
