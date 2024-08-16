%%
clear;
close all;
clc;
run('../../src/Set_Paths');

suppress_figures = true;

M1 = load('fem_matrices.mat', 'mass_matrix').mass_matrix;
M = kron(eye(2), M1);
S1 = load('fem_matrices.mat', 'stiffness_matrix').stiffness_matrix;
S = kron(eye(2), S1);
load('OptimizationSolution.mat');

data_interface = MD_Data_Interface_Transient_ADR_2D();
data_interface.Load_Data();

t = load('OpInf_Training_Data.mat', 't').t;
T = t(end);
n_t = length(t);
n_y = solver.n_y;
n_q = solver.n_q;
beta_t = 50; % Need to look at this more closely
beta_i = 1.e5; % Need to look at this more closely
transient_prior_cov = MD_Transient_Prior_Covariance_Sabl(beta_t, beta_i, T, n_t, n_y);

alpha_u = (1 / 1)^2; % Need to look at this more closely
u_prior_interface = MD_Transient_Elliptic_u_Prior_Interface_Transient_ADR_2D(alpha_u, transient_prior_cov, M, S);

beta_t = 50; % Need to look at this more closely
beta_i = 1.e5; % Need to look at this more closely
transient_prior_cov_z = MD_Transient_Prior_Covariance_Sabl(beta_t, beta_i, t(end - 1), n_t - 1, n_q);
z_prior_interface = MD_z_Prior_Interface_Transient_ADR_2D(transient_prior_cov_z, n_q);

%%
num_prior_samples = 100;
md_prior_sampling = MD_Prior_Sampling(data_interface, u_prior_interface, z_prior_interface);

delta_samples = md_prior_sampling.Prior_Discrepancy_Samples_at_z_opt(num_prior_samples);
if ~suppress_figures

    k = 1;
    u = reshape(delta_samples(:, k), n_y, n_t);
    solver.Animate_Solution(u);

end

%%
z_lofi = sqrt(Q_rom(:));
z = zeros(length(z_lofi), 1);
z(:, 1) = z_lofi + 1;
if ~suppress_figures
    for j = 1:size(z, 2)
        Q = reshape(z(:, j), n_q, n_t - 1).^2;
        figure;
        semilogy(t(2:end), Q);
        title('Optimal controls');
    end
end

delta_prior_samples = md_prior_sampling.Prior_Discrepancy_Samples(z, num_prior_samples);
if ~suppress_figures

    k = 1;
    u = reshape(delta_prior_samples{k}, n_y, n_t);
    solver.Animate_Solution(u);

end

%%
md_post_sampling = MD_Posterior_Sampling(data_interface, u_prior_interface, z_prior_interface);

alpha_d = 1.e-7; % Need to look at this more closely
num_post_samples = 5;
md_post_sampling.Compute_Posterior_Data(alpha_d, num_post_samples);
Z_test = zeros(length(z_lofi), 2);
Z_test(:, 1) = z_lofi;
Z_test(:, 2) = z_lofi + 1.0;
[delta_mean, delta_samples] = md_post_sampling.Posterior_Discrepancy_Samples(Z_test);

if ~suppress_figures

    u = reshape(md_post_sampling.post_data.D(:, 1), n_y, n_t);
    solver.Animate_Solution(u);
    u = reshape(delta_mean{1}, n_y, n_t);
    solver.Animate_Solution(u);

    error = norm(md_post_sampling.post_data.D(:, 1) - delta_mean{1}) / norm(md_post_sampling.post_data.D(:, 1));

    u = reshape(delta_mean{2}, n_y, n_t);
    solver.Animate_Solution(u);
    for k = 1:num_post_samples
        u = reshape(delta_samples{2}(:, k), n_y, n_t);
        solver.Animate_Solution(u);
    end

end

%%

% Need to update MD_Opt_Prob_Interface with a custom ROM class
opt_prob_interface = MD_Opt_Prob_Interface_Transient_ADR_2D(opt, obj_hifi, basis1, basis2, z_lofi);

md_hessian_analysis = MD_Hessian_Analysis(opt_prob_interface, z_prior_interface);
num_evals = 26; % Need to look at this more closely
oversampling = 20;
md_hessian_analysis.Compute_Hessian_GEVP(data_interface.z_opt, num_evals, oversampling);
md_update = MD_Update(md_post_sampling, md_hessian_analysis);

[z_update_mean, z_update_samples] = md_update.Posterior_Update_Samples();

if ~suppress_figures
    Q = reshape(z_lofi, n_q, n_t - 1).^2;
    figure;
    semilogy(t(2:end), Q);
    title('LoFi Optimal controls');

    Q = reshape(z_update_mean, n_q, n_t - 1).^2;
    figure;
    semilogy(t(2:end), Q);
    title('Mean Update Optimal controls');

    Q = reshape(z_hifi, n_q, n_t - 1).^2;
    figure;
    semilogy(t(2:end), Q);
    title('HiFi Optimal controls');
end
