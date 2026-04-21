%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Set up
clear;
close all;
clc;
addpath(genpath('../../src'));
load Optimization_Results.mat;

suppress_figures = true;

n_z = num_space_control_nodes * (n_t - 1);
obj = Adv_Diff_Gaussian_Source_Objective(n_y, n_z, T, n_t, num_space_control_nodes, reg_coeff);
con_lofi = Adv_Diff_Gaussian_Source_Constraint(n_y, n_z, T, n_t, num_space_control_nodes, diff_coeff, vel_coeff_lofi);
opt_lofi = Reduced_Space_Optimization(obj, con_lofi);
con_hifi = Adv_Diff_Gaussian_Source_Constraint(n_y, n_z, T, n_t, num_space_control_nodes, diff_coeff, vel_coeff_hifi);
opt_hifi = Reduced_Space_Optimization(obj, con_hifi);
x = con_lofi.x;

data_interface = MD_Data_Interface_Adv_Diff();
data_interface.Load_Data();

hyperparams = MD_Hyperparameters_Transient_Adv_Diff(data_interface, n_y);

hyperparams.beta_t = 50;
transient_prior_cov = MD_Transient_Prior_Covariance_Sabl(hyperparams, T, n_t, n_y);

alpha_u = (1 / 1)^2;
spatial_u_prior_interface = MD_Elliptic_u_Prior_Interface_Adv_Diff(alpha_u, opt_lofi);
u_prior_interface = MD_Transient_Elliptic_u_Prior_Interface(spatial_u_prior_interface, transient_prior_cov);
z_prior_interface = MD_z_Prior_Interface_Adv_Diff(obj);

%%
num_prior_samples = 100;
md_prior_sampling = MD_Prior_Sampling(data_interface, u_prior_interface, z_prior_interface);

delta_samples = md_prior_sampling.Prior_Discrepancy_Samples_at_z_opt(num_prior_samples);
if ~suppress_figures

    k = 1;
    Plot_State(x, delta_samples(:, k));

    terminal_delta = delta_samples((n_y * (n_t - 1) + 1):end, :);
    figure;
    plot(x, terminal_delta, 'LineWidth', 3, 'color', [.9, .9, .9]);
    plot(x, terminal_delta(:, 1:10), 'LineWidth', 3);
end

%%
z = zeros(length(z_lofi), 1);
z(:, 1) = z_lofi + 1;
if ~suppress_figures
    for j = 1:size(z, 2)
        Plot_Control(x, z(:, j), con_lofi);
    end
end

delta_prior_samples = md_prior_sampling.Prior_Discrepancy_Samples(z, num_prior_samples);
if ~suppress_figures
    for k = 1:10
        Plot_State(x, delta_prior_samples{k});
    end
end

%%
md_post_sampling = MD_Posterior_Sampling(data_interface, u_prior_interface, z_prior_interface);

alpha_d = 1.e-7;
num_post_samples = 100;
md_post_sampling.Compute_Posterior_Data(alpha_d, num_post_samples);
Z_test = zeros(length(z_lofi), 2);
Z_test(:, 1) = z_lofi;
Z_test(:, 2) = z_lofi + 1.0;
[delta_mean, delta_samples] = md_post_sampling.Posterior_Discrepancy_Samples(Z_test);

if ~suppress_figures

    Plot_State(x, md_post_sampling.post_data.D(:, 1));
    Plot_State(x, delta_mean{1});
    error = State_Norm(md_post_sampling.post_data.D(:, 1) - delta_mean{1}, con_lofi) / State_Norm(md_post_sampling.post_data.D(:, 1), con_lofi);

    Plot_State(x, delta_mean{2});
    for k = 1:num_post_samples
        Plot_State(x, delta_samples{2}(:, k));
    end

end

%%
opt_prob_interface = MD_Opt_Prob_Interface_Sabl(opt_lofi, data_interface);
md_hessian_analysis = MD_Hessian_Analysis(opt_prob_interface, z_prior_interface);
num_evals = 26;
oversampling = 20;
md_hessian_analysis.Compute_Hessian_GEVP(data_interface.z_opt, num_evals, oversampling);
md_update = MD_Update(md_post_sampling, md_hessian_analysis);

[z_update_mean, z_update_samples] = md_update.Posterior_Update_Samples();

if ~suppress_figures
    Plot_Control(x, z_lofi, con_lofi);
    Plot_Control(x, z_update_mean, con_lofi);
    Plot_Control(x, z_hifi, con_lofi);
end

%%
u_update = con_hifi.State_Solve(z_update_mean);
u_lofi = con_hifi.State_Solve(z_lofi);
lofi_state_error = State_Norm(u_lofi - u_hifi, con_lofi) / State_Norm(u_hifi, con_lofi);
update_state_error = State_Norm(u_update - u_hifi, con_lofi) / State_Norm(u_hifi, con_lofi);
