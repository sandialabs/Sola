%%
clear;
close all;
clc;
run('../../src/Set_Paths');
rng(13424);

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
u_hyperparam_interface_cell = cell(2, 1);
u_spatial_prior_interface_cell = cell(2, 1);
u_transient_prior_cov = cell(2, 1);
u_prior_interface_cell = cell(2, 1);

k = 1;
u_hyperparam_interface_cell{k} = MD_u_Hyperparameter_Interface_Transient_ADR_2D(true, false, true, k, solver, t);
u_hyperparam_interface_cell{k}.gsvd_num_sing_vals = 1000;
u_hyperparam_interface_cell{k}.gsvd_oversampling = 20;
u_hyperparam_interface_cell{k}.time_variance_inflation = .05;
u_hyperparam_interface_cell{k}.beta_u = .1;
u_spatial_prior_interface_cell{k} = MD_Numeric_Laplacian_u_Prior_Interface(S, M, data_interface, u_hyperparam_interface_cell{k});
u_transient_prior_cov{k} = MD_Transient_Prior_Covariance_Sabl(data_interface, u_hyperparam_interface_cell{k}, T, n_t, n_y / 2);
u_prior_interface_cell{k} = MD_Transient_Elliptic_u_Prior_Interface(data_interface, u_spatial_prior_interface_cell{k}, u_transient_prior_cov{k});

k = 2;
u_hyperparam_interface_cell{k} = MD_u_Hyperparameter_Interface_Transient_ADR_2D(true, false, true, k, solver, t);
u_hyperparam_interface_cell{k}.gsvd_num_sing_vals = 1000;
u_hyperparam_interface_cell{k}.gsvd_oversampling = 20;
u_hyperparam_interface_cell{k}.time_variance_inflation = .2;
u_hyperparam_interface_cell{k}.beta_u = .2;
u_spatial_prior_interface_cell{k} = MD_Numeric_Laplacian_u_Prior_Interface(S, M, data_interface, u_hyperparam_interface_cell{k});
u_transient_prior_cov{k} = MD_Transient_Prior_Covariance_Sabl(data_interface, u_hyperparam_interface_cell{k}, T, n_t, n_y / 2);
u_prior_interface_cell{k} = MD_Transient_Elliptic_u_Prior_Interface(data_interface, u_spatial_prior_interface_cell{k}, u_transient_prior_cov{k});

u_prior_interface = MD_Multi_State_u_Prior_Interface(data_interface, u_prior_interface_cell, u_hyperparam_interface_cell);

%%
num_state_solves = 50;
z_hyperparam_interface = MD_z_Hyperparameter_Interface_Transient_ADR_2D(num_state_solves, opt, basis1, basis2);
z_hyperparam_interface.alpha_z = .009;

h = t(2) - t(1);
m = n_t - 1;
M_t = diag(4 * ones(1, m)) + diag(ones(1, m - 1), 1) + diag(ones(1, m - 1), -1);
M_t(1, 1) = .5 * M_t(1, 1);
M_t(end, end) = .5 * M_t(end, end);
M_t = (1 / 6) * h * M_t;
S_t = diag(2 * ones(1, m)) + (-1) * diag(ones(1, m - 1), 1) + (-1) * diag(ones(1, m - 1), -1);
S_t(1, 1) = .5 * S_t(1, 1);
S_t(end, end) = .5 * S_t(end, end);
S_t = (1 / h) * S_t;
z_prior_interface = MD_Transient_Vector_z_Prior_Interface(S_t, M_t, n_q, data_interface, z_hyperparam_interface, u_prior_interface);

%%
% Uncomment the code below to test hyper-parameter values
% num_prior_samples = 50;
% md_prior_sampling = MD_Prior_Sampling(data_interface, u_prior_interface, z_prior_interface);

%%
% md_prior_sampling.Generate_Prior_Discrepancy_z_opt_Sample_Data(num_prior_samples);
% md_prior_vis = MD_Prior_Visualization(md_prior_sampling);
% md_prior_vis.Visualization_for_Prior_Time_Evolution(1,false);
% md_prior_vis.Visualization_for_Prior_Time_Evolution(2,false);
% md_prior_vis.Visualization_for_Prior_Discrepancy_at_z_opt(1);
% md_prior_vis.Visualization_for_Prior_Discrepancy_at_z_opt(2);

% md_prior_vis.prior_sampling.z_pert_subsample_factor = 20;
% md_prior_vis.prior_sampling.Generate_Prior_Discrepancy_z_pert_Sample_Data();
% md_prior_vis.Visualization_for_Prior_Discrepancy_at_z_pert(1);
% md_prior_vis.Visualization_for_Prior_Discrepancy_at_z_pert(2);

%%
md_post_sampling = MD_Posterior_Sampling(data_interface, u_prior_interface, z_prior_interface);

alpha_d = (1.e-2) * mean([u_prior_interface.u_prior_interface_cell{1}.u_hyperparam_interface.alpha_d, u_prior_interface.u_prior_interface_cell{2}.u_hyperparam_interface.alpha_d]);
num_post_samples = 50;
md_post_sampling.Compute_Posterior_Data(alpha_d, num_post_samples);

z_lofi = data_interface.Load_Optimal_z();
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

obj_update_samples = zeros(num_post_samples, 1);
for k = 1:num_post_samples
    disp(['Computer posterior sample ', num2str(k)]);
    Q_update_sample = reshape(z_update_samples(:, k), n_q, n_t - 1);
    pp = pchip(t, [Q_update_sample(:, 1), Q_update_sample].^2);
    controller = @(tt) ppval(pp, tt);
    Y_update_sample = solver.State_Solve(controller, t).NodalSolution;
    u_tmp1 = reshape(Y_update_sample(:, 1, :), [], n_t);
    u_tmp2 = reshape(Y_update_sample(:, 2, :), [], n_t);
    utmp = [u_tmp1; u_tmp2];
    u_update_sample = utmp(:);
    obj_update_samples(k) = obj_hifi.J(u_update_sample, z_update_samples(:, k));
end

%%
% save('MD_Results.mat', 'z_lofi', 'z_update_mean', 'z_update_samples', 'Y_lofi', 'Y_update_mean','obj_update_samples','obj_update_mean','obj_lofi');
