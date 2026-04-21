%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
close all;
clc;
addpath(genpath('../../src'));
rng(12342);

m = load('Optimization_Results.mat', 'm').m;
vel_coeff = load('Optimization_Results.mat', 'vel_coeff').vel_coeff;
reg_coeff = load('Optimization_Results.mat', 'reg_coeff').reg_coeff;
diff_coeff = load('Optimization_Results.mat', 'diff_coeff').diff_coeff;
n_r = size(diff_coeff, 2);

obj = Adv_Diff_Objective(m, reg_coeff);
cons_hifi = cell(n_r, 1);
cons_lofi = cell(n_r, 1);
for k = 1:n_r
    cons_hifi{k} = Adv_Diff(m, vel_coeff, diff_coeff(k));
    cons_lofi{k} = Diff(cons_hifi{k});
end
opt_lofi = Reduced_Space_Optimization_Under_Uncertainty(obj, cons_lofi);

md_ouu_data_interface = MD_Data_Interface_Diff();
u_hyperparam_interface = load('Hyperparameter_Interfaces.mat', 'u_hyperparam_interface').u_hyperparam_interface;
z_hyperparam_interface = load('Hyperparameter_Interfaces.mat', 'z_hyperparam_interface').z_hyperparam_interface;

x = cons_lofi{1}.x;
S = cons_lofi{1}.S;
M = cons_lofi{1}.M;

us_prior_interface = MD_Numeric_Laplacian_u_Prior_Interface(S, M, md_ouu_data_interface, u_hyperparam_interface);
ensemble_weighting = MD_OUU_Ensemble_Weighting_Matrix(md_ouu_data_interface, us_prior_interface);
u_prior_interface = MD_OUU_u_Prior_Interface(us_prior_interface, md_ouu_data_interface, ensemble_weighting);
z_prior_interface = MD_Numeric_Laplacian_z_Prior_Interface(S, M, md_ouu_data_interface, z_hyperparam_interface, u_prior_interface);

%%
md_post_sampling = MD_Posterior_Sampling(md_ouu_data_interface, u_prior_interface, z_prior_interface);
alpha_d = u_hyperparam_interface.alpha_d;
num_post_samples = 100;
md_post_sampling.Compute_Posterior_Data(alpha_d, num_post_samples);
Z_test = zeros(m, 3);
Z_test(:, 1:2) = md_ouu_data_interface.Z;
Z_test(:, 3) = 0.5 * md_ouu_data_interface.Z(:, 1) + 0.5 * md_ouu_data_interface.Z(:, 2) + 0.7 * x;
[delta_mean, delta_samples] = md_post_sampling.Posterior_Discrepancy_Samples(Z_test);

%%
delta_mean_z1 = md_ouu_data_interface.Reshape_State_to_Mat(delta_mean{1});
n_s = size(delta_samples{1}, 2);
delta_samples_z1 = zeros(m, n_s, n_r);
for k = 1:n_s
    delta_samples_z1(:, k, :) = md_ouu_data_interface.Reshape_State_to_Mat(delta_samples{1}(:, k));
end
discrepancy_z1 = zeros(m, n_r);
for j = 1:n_r
    discrepancy_z1(:, j) = cons_hifi{j}.State_Solve(Z_test(:, 1)) - cons_lofi{j}.State_Solve(Z_test(:, 1));
end

delta_mean_z3 = md_ouu_data_interface.Reshape_State_to_Mat(delta_mean{3});
delta_samples_z3 = zeros(m, n_s, n_r);
for k = 1:n_s
    delta_samples_z3(:, k, :) = md_ouu_data_interface.Reshape_State_to_Mat(delta_samples{3}(:, k));
end
discrepancy_z3 = zeros(m, n_r);
for j = 1:n_r
    discrepancy_z3(:, j) = cons_hifi{j}.State_Solve(Z_test(:, 3)) - cons_lofi{j}.State_Solve(Z_test(:, 3));
end

%%
r = 15;
xi1 = 1;
xi2 = 2;
xi3 = 3;

figure;
hold on;
plot(x, discrepancy_z1(:, xi1) - discrepancy_z1(:, xi2), "Color", 'magenta', 'LineWidth', 3);
plot(x, delta_samples_z1(:, :, xi1) - delta_samples_z1(:, :, xi2), "Color", [.9, .9, .9], 'LineWidth', 3);
plot(x, discrepancy_z1(:, xi1) - discrepancy_z1(:, xi2), "Color", 'magenta', 'LineWidth', 3);
legend({'$(S(z_1,\xi_1)-\tilde{S}(z_1,\xi_1))-(S(z_1,\xi_2)-\tilde{S}(z_1,\xi_2))$', '$\delta(z_1,\theta_1)-\delta(z_1,\theta_2)$'}, 'Interpreter', 'latex');
ylim([-r, r]);
set(gca, 'fontsize', 20);
% exportgraphics(gcf, 'z_1_xi_12_corr.eps', 'BackgroundColor', 'none', 'ContentType', 'vector');

figure;
hold on;
plot(x, discrepancy_z1(:, xi1) - discrepancy_z1(:, xi3), "Color", 'magenta', 'LineWidth', 3);
plot(x, delta_samples_z1(:, :, xi1) - delta_samples_z1(:, :, xi3), "Color", [.9, .9, .9], 'LineWidth', 3);
plot(x, discrepancy_z1(:, xi1) - discrepancy_z1(:, xi3), "Color", 'magenta', 'LineWidth', 3);
legend({'$(S(z_1,\xi_3)-\tilde{S}(z_1,\xi_3))-(S(z_1,\xi_3)-\tilde{S}(z_1,\xi_3))$', '$\delta(z_1,\theta_1)-\delta(z_1,\theta_3)$'}, 'Interpreter', 'latex');
ylim([-r, r]);
set(gca, 'fontsize', 20);
% exportgraphics(gcf, 'z_1_xi_13_corr.eps', 'BackgroundColor', 'none', 'ContentType', 'vector');

figure;
hold on;
plot(x, discrepancy_z3(:, xi1) - discrepancy_z3(:, xi2), "Color", 'magenta', 'LineWidth', 3);
plot(x, delta_samples_z3(:, :, xi1) - delta_samples_z3(:, :, xi2), "Color", [.9, .9, .9], 'LineWidth', 3);
plot(x, discrepancy_z3(:, xi1) - discrepancy_z3(:, xi2), "Color", 'magenta', 'LineWidth', 3);
legend({'$(S(z_3,\xi_1)-\tilde{S}(z_3,\xi_1))-(S(z_3,\xi_2)-\tilde{S}(z_3,\xi_2))$', '$\delta(z_3,\theta_1)-\delta(z_3,\theta_2)$'}, 'Interpreter', 'latex');
ylim([-r, r]);
set(gca, 'fontsize', 20);
% exportgraphics(gcf, 'z_3_xi_12_corr.eps', 'BackgroundColor', 'none', 'ContentType', 'vector');

figure;
hold on;
plot(x, discrepancy_z3(:, xi1) - discrepancy_z3(:, xi3), "Color", 'magenta', 'LineWidth', 3);
plot(x, delta_samples_z3(:, :, xi1) - delta_samples_z3(:, :, xi3), "Color", [.9, .9, .9], 'LineWidth', 3);
plot(x, discrepancy_z3(:, xi1) - discrepancy_z3(:, xi3), "Color", 'magenta', 'LineWidth', 3);
legend({'$(S(z_3,\xi_1)-\tilde{S}(z_3,\xi_1))-(S(z_3,\xi_3)-\tilde{S}(z_3,\xi_3))$', '$\delta(z_3,\theta_1)-\delta(z_3,\theta_3)$'}, 'Interpreter', 'latex');
ylim([-r, r]);
set(gca, 'fontsize', 20);
% exportgraphics(gcf, 'z_3_xi_13_corr.eps', 'BackgroundColor', 'none', 'ContentType', 'vector');
