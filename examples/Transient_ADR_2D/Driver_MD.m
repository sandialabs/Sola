%%
clear;
close all;
clc;
run('../../src/Set_Paths');

suppress_figures = false;

M = load('fem_matrices.mat', 'mass_matrix').mass_matrix;
S = load('fem_matrices.mat', 'stiffness_matrix').stiffness_matrix;

load('OptimizationSolution.mat');
basis1.Set_Reduced_Dimension_From_Residual_Energy(residual_energies(1));
basis2.Set_Reduced_Dimension_From_Residual_Energy(residual_energies(1));

data_interface = MD_Data_Interface_Transient_ADR_2D(solver);

t = load('OpInf_Training_Data.mat', 't').t;
T = t(end);
n_t = length(t);
n_y = solver.n_y;
n_q = solver.n_q;

%%
ui_hyperparams = cell(2,1);
ui_spatial_prior_interface = cell(2,1);
ui_transient_prior_cov = cell(2,1);
ui_prior_interface = cell(2,1);
for k = 1:2
    ui_hyperparams{k} = MD_u_Hyperparameters_Transient_ADR_2D(data_interface, solver, t, k);
    ui_hyperparams{k}.gsvd_num_sing_vals = 1000;
    ui_hyperparams{k}.gsvd_oversampling = 20;
    ui_spatial_prior_interface{k} = MD_Numeric_Laplacian_u_Prior_Interface(S,M,ui_hyperparams{k});
    ui_transient_prior_cov{k} = MD_Transient_Prior_Covariance_Sabl(ui_hyperparams{k}, T, n_t, n_y/2);
    ui_prior_interface{k} = MD_Transient_Elliptic_u_Prior_Interface(ui_spatial_prior_interface{k}, ui_transient_prior_cov{k});
end
u_prior_interface = MD_Multi_State_u_Prior_Interface(ui_prior_interface);

%%
z_hyperparams = struct;
z_hyperparams.beta_t = 1.0;
transient_prior_cov_z = MD_Transient_Prior_Covariance_Sabl(z_hyperparams, t(end - 1), n_t - 1, n_q);
z_prior_interface = MD_z_Prior_Interface_Transient_ADR_2D(transient_prior_cov_z, n_q);

time_series_samples_z = transient_prior_cov_z.Sample_Time_Series(10);
time_series_samples_u = ui_transient_prior_cov.Sample_Time_Series(10);

%%
num_prior_samples = 10;
md_prior_sampling = MD_Prior_Sampling(data_interface, u_prior_interface, z_prior_interface);

prior_delta_samples_z_opt = md_prior_sampling.Prior_Discrepancy_Samples_at_z_opt(num_prior_samples);
if ~suppress_figures

    u = data_interface.D;
    u_tmp = reshape(u, n_y, n_t);
    u1 = u_tmp(1:(n_y / 2), :);
    u2 = u_tmp((n_y / 2 + 1):end, :);
    data1_min = min(u1, [], 1);
    data1_max = max(u1, [], 1);
    data2_min = min(u2, [], 1);
    data2_max = max(u2, [], 1);

    delta1_min = zeros(n_t, num_prior_samples);
    delta1_max = zeros(n_t, num_prior_samples);
    delta2_min = zeros(n_t, num_prior_samples);
    delta2_max = zeros(n_t, num_prior_samples);
    for k = 1:num_prior_samples
        u = prior_delta_samples_z_opt(:, k);
        u_tmp = reshape(u, n_y, n_t);
        u1 = u_tmp(1:(n_y / 2), :);
        u2 = u_tmp((n_y / 2 + 1):end, :);
        delta1_min(:, k) = min(u1, [], 1);
        delta1_max(:, k) = max(u1, [], 1);
        delta2_min(:, k) = min(u2, [], 1);
        delta2_max(:, k) = max(u2, [], 1);
    end

    figure;
    hold on;
    for k = 1:num_prior_samples
        plot(t, delta1_min(:, k), 'LineWidth', 3, 'color', .9 * [1, 1, 1]);
        plot(t, delta1_max(:, k), 'LineWidth', 3, 'color', .9 * [1, 1, 1]);
    end
    plot(t, data1_min, 'LineWidth', 3, 'color', 'black');
    plot(t, data1_max, 'LineWidth', 3, 'color', 'black');
    title('Discrepancy 1 Magnitude');

    figure;
    hold on;
    for k = 1:num_prior_samples
        plot(t, delta2_min(:, k), 'LineWidth', 3, 'color', .9 * [1, 1, 1]);
        plot(t, delta2_max(:, k), 'LineWidth', 3, 'color', .9 * [1, 1, 1]);
    end
    plot(t, data2_min, 'LineWidth', 3, 'color', 'black');
    plot(t, data2_max, 'LineWidth', 3, 'color', 'black');
    title('Discrepancy 2 Magnitude');

    k = 1;
    time_step = 1;

    u = reshape(data_interface.D, n_y, n_t);
    figure;
    pdeplot(solver.model.Mesh, XYData = u(1:(n_y / 2), time_step), colormap = 'parula');
    u = reshape(prior_delta_samples_z_opt(:, k), n_y, n_t);
    figure;
    pdeplot(solver.model.Mesh, XYData = u(1:(n_y / 2), time_step), colormap = 'parula');

    u = reshape(data_interface.D, n_y, n_t);
    figure;
    pdeplot(solver.model.Mesh, XYData = u((n_y / 2 + 1):end, time_step), colormap = 'parula');
    u = reshape(prior_delta_samples_z_opt(:, k), n_y, n_t);
    figure;
    pdeplot(solver.model.Mesh, XYData = u((n_y / 2 + 1):end, time_step), colormap = 'parula');

end

%%
z_lofi = Q_rom(:);
z = zeros(length(z_lofi), 1);
z(:, 1) = z_lofi .* 1.2;
if ~suppress_figures
    for j = 1:size(z, 2)
        Q = reshape(z(:, j), n_q, n_t - 1).^2;
        figure;
        semilogy(t(2:end), Q);
        title('Optimal controls');
    end
end

prior_delta_samples = md_prior_sampling.Prior_Discrepancy_Samples(z, num_prior_samples);
if ~suppress_figures

    k = 2;
    u = reshape(prior_delta_samples{k}, n_y, n_t);
    time_step = 5;
    figure;
    pdeplot(solver.model.Mesh, XYData = u((n_y / 2 + 1):end, time_step), colormap = 'parula');

    u = reshape(prior_delta_samples_z_opt(:, k), n_y, n_t);
    figure;
    pdeplot(solver.model.Mesh, XYData = u((n_y / 2 + 1):end, time_step), colormap = 'parula');

    u = reshape(data_interface.u_opt, n_y, n_t);
    figure;
    pdeplot(solver.model.Mesh, XYData = u((n_y / 2 + 1):end, time_step), colormap = 'parula');

end

%%
md_post_sampling = MD_Posterior_Sampling(data_interface, u_prior_interface, z_prior_interface);

alpha_d = 1.e-8;
num_post_samples = 1;
md_post_sampling.Compute_Posterior_Data(alpha_d, num_post_samples);
Z_test = zeros(length(z_lofi), 2);
Z_test(:, 1) = z_lofi;
Z_test(:, 2) = z_lofi .* 1.2;
[post_delta_mean, post_delta_samples] = md_post_sampling.Posterior_Discrepancy_Samples(Z_test);

data_fit_error_mean = norm(md_post_sampling.post_data.D(:, 1) - post_delta_mean{1}) / norm(md_post_sampling.post_data.D(:, 1));

data_fit_error_samples_1 = zeros(num_post_samples, 1);
data_fit_error_samples_2 = zeros(num_post_samples, 1);
for k = 1:num_post_samples
    data_fit_error_samples_1(k) = norm(md_post_sampling.post_data.D(:, 1) - post_delta_samples{1}(:, k)) / norm(md_post_sampling.post_data.D(:, 1));
    data_fit_error_samples_2(k) = norm(md_post_sampling.post_data.D(:, 1) - post_delta_samples{2}(:, k)) / norm(md_post_sampling.post_data.D(:, 1));
end

%%
opt_prob_interface = MD_Opt_Prob_Interface_Transient_ADR_2D(opt, obj_hifi, basis1, basis2, z_lofi);

md_hessian_analysis = MD_Hessian_Analysis(opt_prob_interface, z_prior_interface);
num_evals = 15;
oversampling = 20;
md_hessian_analysis.Compute_Hessian_GEVP(data_interface.z_opt, num_evals, oversampling);
md_update = MD_Update(md_post_sampling, md_hessian_analysis);

[z_update_mean, z_update_samples] = md_update.Posterior_Update_Samples();

if ~suppress_figures
    Q_lofi = reshape(z_lofi, n_q, n_t - 1);
    figure;
    semilogy(t(2:end), abs(Q_lofi));
    title('LoFi Optimal controls');

    Q_update_mean = reshape(z_update_mean, n_q, n_t - 1);
    figure;
    semilogy(t(2:end), abs(Q_update_mean));
    title('Mean Update Optimal controls');
end

%%
Q_lofi = reshape(z_lofi, n_q, n_t - 1);
pp = pchip(t, [Q_lofi(:, 1), Q_lofi].^2);
controller = @(tt) ppval(pp, tt);
Y_lofi = solver.State_Solve(controller, t).NodalSolution;
u_tmp1 = reshape(Y_lofi(:, 1, :), [], n_t);
u_tmp2 = reshape(Y_lofi(:, 2, :), [], n_t);
utmp = [u_tmp1; u_tmp2];
u_lofi = utmp(:);
obj_lofi = obj_hifi.J(u_lofi, z_lofi);

Q_update_mean = reshape(z_update_mean, n_q, n_t - 1);
pp = pchip(t, [Q_update_mean(:, 1), Q_update_mean].^2);
controller = @(tt) ppval(pp, tt);
Y_update_mean = solver.State_Solve(controller, t).NodalSolution;
u_tmp1 = reshape(Y_update_mean(:, 1, :), [], n_t);
u_tmp2 = reshape(Y_update_mean(:, 2, :), [], n_t);
utmp = [u_tmp1; u_tmp2];
u_update_mean = utmp(:);
obj_update_mean = obj_hifi.J(u_update_mean, z_update_mean);

%%
save('MD_Results.mat', 'z_lofi', 'z_update_mean', 'z_update_samples', 'Y_lofi', 'Y_update_mean');
