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
rng(1232453);

working_path = pwd;
write_path = '~/Documents/dasco/papers/Model_Discrepancy_Hyperparameters/figures/';

alpha_u_pert = 100;
beta_u_pert = 100;
alpha_z_pert = 100;
beta_z_pert = 10;

con_hifi = Adv_Diff_Constraint(m, diff_coeff, vel_coeff, robin_coeff);
con_lofi = Diff_Constraint(con_hifi);
x = con_lofi.x;

num_samples = 2;

E_z = (1.e-1) * con_hifi.S + con_hifi.M;

z = linsolve(E_z, sqrtm(con_hifi.M) * randn(m, num_samples - 1));
z_lofi_norm = sqrt(z_lofi' * con_hifi.M * z_lofi);
z_norm = sqrt(diag(z' * con_hifi.M * z));
for k = 1:(num_samples - 1)
    z(:, k) = z_lofi + .3 * z_lofi_norm * z(:, k) / z_norm(k);
end
Z = zeros(m, num_samples);
Z(:, 1) = z_lofi;
Z(:, 2:end) = z;
D = con_hifi.State_Solve(Z) - con_lofi.State_Solve(Z);

data_interface = MD_Data_Interface_Discrepancy_Calibration(z_lofi, u_lofi, Z, D);
data_centering = true;

u_hyperparam_interface = MD_u_Hyperparameter_Interface_Discrepancy_Calibration(x, data_centering);
u_prior_interface = MD_Numeric_Laplacian_u_Prior_Interface(con_hifi.S, con_hifi.M, data_interface, u_hyperparam_interface);

num_state_solves = 100;
z_hyperparam_interface = MD_z_Hyperparameter_Interface_Discrepancy_Calibration(num_state_solves, x, con_lofi);
z_prior_interface = MD_Numeric_Laplacian_z_Prior_Interface(con_hifi.S, con_hifi.M, data_interface, z_hyperparam_interface, u_prior_interface);

alpha_u = alpha_u_pert / 100 * u_hyperparam_interface.alpha_u;
u_prior_interface.Set_alpha_u(alpha_u);
beta_u = beta_u_pert / 100 * u_hyperparam_interface.beta_u;
u_prior_interface.Set_beta_u(beta_u);

alpha_z = alpha_z_pert / 100 * z_hyperparam_interface.alpha_z;
z_prior_interface.Set_alpha_z(alpha_z);
beta_z = beta_z_pert / 100 * z_hyperparam_interface.beta_z;
z_prior_interface.Set_beta_z(beta_z);

%%
md_post_sampling = MD_Posterior_Sampling(data_interface, u_prior_interface, z_prior_interface);
alpha_d = u_hyperparam_interface.alpha_d;
num_post_samples = 100;
md_post_sampling.Compute_Posterior_Data(alpha_d, num_post_samples);

obj = Adv_Diff_Objective(m, reg_coeff);
opt_lofi = Reduced_Space_Optimization(obj, con_lofi);
opt_prob_interface = MD_Opt_Prob_Interface_Sabl(opt_lofi, data_interface);
md_hessian_analysis = MD_Hessian_Analysis(opt_prob_interface, z_prior_interface);
md_update = MD_Update(md_post_sampling, md_hessian_analysis);

[z_update_mean, z_update_samples] = md_update.Posterior_Update_Samples();

figure;
hold on;
plot(x, z_lofi, 'color', 'black', 'LineWidth', 3);
plot(x, z_hifi, 'color', 'cyan', 'LineWidth', 3);
plot(x, z_update_mean, '--', 'color', 'red', 'LineWidth', 3);
for k = 1:num_post_samples
    plot(x, z_update_samples(:, k), 'color', [.9, .9, .9], 'LineWidth', 3);
end
plot(x, md_post_sampling.z_opt, 'color', 'black', 'LineWidth', 3);
plot(x, z_hifi, 'color', 'cyan', 'LineWidth', 3);
plot(x, z_update_mean, '--', 'color', 'red', 'LineWidth', 3);
xlabel('Spatial Input');
ylabel('Optimal Controller');
ylim([.5, 2]);
legend({'Low-fidelity Control', 'High-fidelity Control', 'Posterior Mean', 'Posterior Samples'}, 'Position', [0.1696    0.2036    0.4125    0.2607]);
set(gca, 'fontsize', 24);
cd(write_path);
saveas(gca, ['posterior_z_samples_', num2str(alpha_u_pert), '_', num2str(beta_u_pert), '_', num2str(alpha_z_pert), '_', num2str(beta_z_pert)], 'epsc');
cd(working_path);
