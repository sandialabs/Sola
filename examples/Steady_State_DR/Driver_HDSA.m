%% Set up
clear;
close all;
clc;
addpath(genpath('../../src'));
load Optimization_Results.mat;

obj = Diff_React_Objective(m, reg_coeff);
con_lofi = Diff_React_Constraint(m, diff_coeff, react_coeff);
opt_lofi = Reduced_Space_Optimization(obj, con_lofi);
x = con_lofi.x;

%%
alpha_u = 2^2;
alpha_z = 1.e-10;
md_interface = Diff_React_HDSA(opt_lofi, alpha_u, alpha_z);

%%
num_prior_samples = 500;
md_prior_sampling = HDSA_MD_Prior_Sampling(md_interface);

prior_delta_samples = md_prior_sampling.Prior_Discrepancy_Samples_at_z_opt(num_prior_samples);
figure;
hold on;
plot(x, prior_delta_samples, 'LineWidth', 3, 'color', [.9, .9, .9]);
plot(x, prior_delta_samples(:, 1:10), 'LineWidth', 3);
set(gca, 'fontsize', 18);

z_samples = md_prior_sampling.Prior_z_Samples(10);
figure;
hold on;
plot(x, z_samples, 'LineWidth', 3);
set(gca, 'fontsize', 18);

%%
md_update = HDSA_MD_Update(md_interface);

alpha_d = 1.e-4;
num_post_samples = 500;
md_update.Compute_Posterior_Data(alpha_d, num_post_samples);
Z_test = zeros(m, 3);
Z_test(:, 1:2) = Z;
Z_test(:, 3) = 300 * exp(-10 * (x - 0.5).^2);
[delta_mean, delta_samples] = md_update.Posterior_Discrepancy_Samples(Z_test);

figure;
hold on;
plot(x, Z_test(:, 1), 'LineWidth', 3);
plot(x, Z_test(:, 2), 'LineWidth', 3);
plot(x, Z_test(:, 3), 'LineWidth', 3);
legend({'z_1', 'z_2', 'z_3'});
set(gca, 'fontsize', 18);

figure;
hold on;
plot(x, md_update.post_data.D(:, 1), 'color', 'black', 'LineWidth', 3);
plot(x, delta_mean{1}, '--', 'color', 'red', 'LineWidth', 3);
for k = 1:num_post_samples
    plot(x, delta_samples{1}(:, k), 'color', [.9, .9, .9], 'LineWidth', 3);
end
plot(x, md_update.post_data.D(:, 1), 'color', 'black', 'LineWidth', 3);
plot(x, delta_mean{1}, '--', 'color', 'red', 'LineWidth', 3);
set(gca, 'fontsize', 18);

figure;
hold on;
plot(x, md_update.post_data.D(:, 2), 'color', 'black', 'LineWidth', 3);
plot(x, delta_mean{2}, '--', 'color', 'red', 'LineWidth', 3);
for k = 1:num_post_samples
    plot(x, delta_samples{2}(:, k), 'color', [.9, .9, .9], 'LineWidth', 3);
end
plot(x, md_update.post_data.D(:, 2), 'color', 'black', 'LineWidth', 3);
plot(x, delta_mean{2}, '--', 'color', 'red', 'LineWidth', 3);
set(gca, 'fontsize', 18);

figure;
hold on;
plot(x, delta_mean{3}, '--', 'color', 'red', 'LineWidth', 3);
for k = 1:num_post_samples
    plot(x, delta_samples{3}(:, k), 'color', [.9, .9, .9], 'LineWidth', 3);
end
plot(x, delta_mean{3}, '--', 'color', 'red', 'LineWidth', 3);
set(gca, 'fontsize', 18);

%%
num_evals = 4;
oversampling = 20;
md_update.Compute_Hessian_GEVP(num_evals, oversampling);
[z_update_mean, z_update_samples] = md_update.Posterior_Update_Samples();

figure;
hold on;
plot(x, z_lofi, 'color', 'black', 'LineWidth', 3);
plot(x, z_hifi, 'color', 'cyan', 'LineWidth', 3);
plot(x, z_update_mean, '--', 'color', 'red', 'LineWidth', 3);
for k = 1:num_post_samples
    plot(x, z_update_samples(:, k), 'color', [.9, .9, .9], 'LineWidth', 3);
end
plot(x, md_update.z_opt, 'color', 'black', 'LineWidth', 3);
plot(x, z_hifi, 'color', 'cyan', 'LineWidth', 3);
plot(x, z_update_mean, '--', 'color', 'red', 'LineWidth', 3);
legend({'Low-fidelity control', 'High-fidelity control', 'Update'});
set(gca, 'fontsize', 18);

%%
con_hifi = Diff_React_HiFi_Constraint(con_lofi);
opt_hifi = Reduced_Space_Optimization(obj, con_hifi);

u_true_lofi = con_hifi.State_Solve(z_lofi);
Jhat_lofi = opt_hifi.Jhat(z_lofi);

u_true_hifi = con_hifi.State_Solve(z_hifi);
Jhat_hifi = opt_hifi.Jhat(z_hifi);

u_true_update = con_hifi.State_Solve(z_update_mean);
Jhat_update = opt_hifi.Jhat(z_update_mean);

u_true_update_samples = zeros(m, num_post_samples);
Jhat_update_samples = zeros(num_post_samples, 1);
for k = 1:num_post_samples
    u_true_update_samples(:, k) = con_hifi.State_Solve(z_update_samples(:, k));
    Jhat_update_samples(k) = opt_hifi.Jhat(z_update_samples(:, k));
end

figure;
hold on;
plot(x, u_true_lofi, 'color', 'black', 'LineWidth', 3);
plot(x, u_true_hifi, 'color', 'cyan', 'LineWidth', 3);
plot(x, u_true_update, '--', 'color', 'red', 'LineWidth', 3);
for k = 1:num_post_samples
    plot(x, u_true_update_samples(:, k), 'color', [.9, .9, .9], 'LineWidth', 3);
end
plot(x, u_true_lofi, 'color', 'black', 'LineWidth', 3);
plot(x, u_true_hifi, 'color', 'cyan', 'LineWidth', 3);
plot(x, u_true_update, '--', 'color', 'red', 'LineWidth', 3);
legend({'Low-fidelity controlled state', 'High-fidelity controlled state', 'Update controlled state'});
set(gca, 'fontsize', 18);

figure;
hold on;
histogram(Jhat_update_samples);
plot([Jhat_update, Jhat_update], [0, 40], 'LineWidth', 3);
plot([Jhat_hifi, Jhat_hifi], [0, 40], 'LineWidth', 3);
plot([Jhat_lofi, Jhat_lofi], [0, 40], 'LineWidth', 3);
xlabel('High-fidelity objective function value');
yticks([]);
set(gca, 'fontsize', 18);

%%
rank_range = (1:100)';
oversampling = 20;
z_update_mean_range = cell(length(rank_range), 1);
z_update_samples_range = cell(length(rank_range), 1);
count = 1;
for k = 1:length(rank_range)
    num_evals = rank_range(k);
    md_update.Compute_Hessian_GEVP(num_evals, oversampling);
    [z_update_mean_range{count}, z_update_samples_range{count}] = md_update.Posterior_Update_Samples();
    count = count + 1;
end

save('HDSA_Results.mat');
