%% Set up
clear;
close all;
clc;
addpath(genpath('../../src'));
load Optimization_Results.mat;
rng(2451423);

N = length(z_lofi) + 1;
obj_hifi = Mass_Spring_Objective_HiFi(T, N);
obj_lofi = Mass_Spring_Objective_LoFi(T, N);
con_hifi = Mass_Spring_Coupled(T, N);
con_lofi = Mass_Spring_LoFi(con_hifi);
opt_hifi = Reduced_Space_Optimization(obj_hifi, con_hifi);
opt_lofi = Reduced_Space_Optimization(obj_lofi, con_lofi);
t = con_hifi.t_mesh;

%%
alpha_u = 1.e4;
alpha_z = 1.e-10;
md_interface = Mass_Spring_HDSA(opt_lofi, alpha_u, alpha_z);

%%
num_prior_samples = 500;
md_prior_sampling = HDSA_MD_Prior_Sampling(md_interface);

prior_delta_samples = md_prior_sampling.Prior_Discrepancy_Samples_at_z_opt(num_prior_samples);
figure;
hold on;
plot(t, prior_delta_samples(1:2:end, :), 'LineWidth', 3, 'color', [.9, .9, .9]);
plot(t, prior_delta_samples(1:2:end, 1:5), 'LineWidth', 3);
xlabel('Time');
ylabel('$x_1$', 'Interpreter', 'latex');
set(gca, 'fontsize', 18);

figure;
hold on;
plot(t, prior_delta_samples(2:2:end, :), 'LineWidth', 3, 'color', [.9, .9, .9]);
plot(t, prior_delta_samples(2:2:end, 1:5), 'LineWidth', 3);
xlabel('Time');
ylabel('$v_1$', 'Interpreter', 'latex');
set(gca, 'fontsize', 18);

z_samples = md_prior_sampling.Prior_z_Samples(5);
figure;
hold on;
plot(t(2:end), z_samples, 'LineWidth', 3);
set(gca, 'fontsize', 18);

%%
md_update = HDSA_MD_Update(md_interface);

alpha_d = 1.e-1;
num_post_samples = 500;
md_update.Compute_Posterior_Data(alpha_d, num_post_samples);
Z_test = zeros(N - 1, 3);
Z_test(:, 1:2) = Z;
Z_test(:, 3) = 500 * ones(N - 1, 1);
[delta_mean, delta_samples] = md_update.Posterior_Discrepancy_Samples(Z_test);

figure;
hold on;
plot(t(2:end), Z_test(:, 1), 'LineWidth', 3);
plot(t(2:end), Z_test(:, 2), 'LineWidth', 3);
plot(t(2:end), Z_test(:, 3), 'LineWidth', 3);
legend({'z_1', 'z_2', 'z_3'});
set(gca, 'fontsize', 18);

figure;
hold on;
plot(t, md_update.post_data.D(1:2:end, 1), 'color', 'black', 'LineWidth', 3);
plot(t, delta_mean{1}(1:2:end), '--', 'color', 'red', 'LineWidth', 3);
for k = 1:num_post_samples
    plot(t, delta_samples{1}(1:2:end, k), 'color', [.9, .9, .9], 'LineWidth', 3);
end
plot(t, md_update.post_data.D(1:2:end, 1), 'color', 'black', 'LineWidth', 3);
plot(t, delta_mean{1}(1:2:end), '--', 'color', 'red', 'LineWidth', 3);
xlabel('Time');
ylabel('$x_1$', 'Interpreter', 'latex');
set(gca, 'fontsize', 18);

figure;
hold on;
plot(t, md_update.post_data.D(2:2:end, 1), 'color', 'black', 'LineWidth', 3);
plot(t, delta_mean{1}(2:2:end), '--', 'color', 'red', 'LineWidth', 3);
for k = 1:num_post_samples
    plot(t, delta_samples{1}(2:2:end, k), 'color', [.9, .9, .9], 'LineWidth', 3);
end
plot(t, md_update.post_data.D(2:2:end, 1), 'color', 'black', 'LineWidth', 3);
plot(t, delta_mean{1}(2:2:end), '--', 'color', 'red', 'LineWidth', 3);
xlabel('Time');
ylabel('$v_1$', 'Interpreter', 'latex');
set(gca, 'fontsize', 18);

figure;
hold on;
plot(t, md_update.post_data.D(1:2:end, 2), 'color', 'black', 'LineWidth', 3);
plot(t, delta_mean{2}(1:2:end), '--', 'color', 'red', 'LineWidth', 3);
for k = 1:num_post_samples
    plot(t, delta_samples{2}(1:2:end, k), 'color', [.9, .9, .9], 'LineWidth', 3);
end
plot(t, md_update.post_data.D(1:2:end, 2), 'color', 'black', 'LineWidth', 3);
plot(t, delta_mean{2}(1:2:end), '--', 'color', 'red', 'LineWidth', 3);
xlabel('Time');
ylabel('$x_1$', 'Interpreter', 'latex');
set(gca, 'fontsize', 18);

figure;
hold on;
plot(t, md_update.post_data.D(2:2:end, 2), 'color', 'black', 'LineWidth', 3);
plot(t, delta_mean{2}(2:2:end), '--', 'color', 'red', 'LineWidth', 3);
for k = 1:num_post_samples
    plot(t, delta_samples{2}(2:2:end, k), 'color', [.9, .9, .9], 'LineWidth', 3);
end
plot(t, md_update.post_data.D(2:2:end, 2), 'color', 'black', 'LineWidth', 3);
plot(t, delta_mean{2}(2:2:end), '--', 'color', 'red', 'LineWidth', 3);
xlabel('Time');
ylabel('$v_1$', 'Interpreter', 'latex');
set(gca, 'fontsize', 18);

figure;
hold on;
plot(t, delta_mean{3}(1:2:end), '--', 'color', 'red', 'LineWidth', 3);
for k = 1:num_post_samples
    plot(t, delta_samples{3}(1:2:end, k), 'color', [.9, .9, .9], 'LineWidth', 3);
end
plot(t, delta_mean{3}(1:2:end), '--', 'color', 'red', 'LineWidth', 3);
xlabel('Time');
ylabel('$x_1$', 'Interpreter', 'latex');
set(gca, 'fontsize', 18);

figure;
hold on;
plot(t, delta_mean{3}(2:2:end), '--', 'color', 'red', 'LineWidth', 3);
for k = 1:num_post_samples
    plot(t, delta_samples{3}(2:2:end, k), 'color', [.9, .9, .9], 'LineWidth', 3);
end
plot(t, delta_mean{3}(2:2:end), '--', 'color', 'red', 'LineWidth', 3);
xlabel('Time');
ylabel('$v_1$', 'Interpreter', 'latex');
set(gca, 'fontsize', 18);

%%
num_evals = 17;
oversampling = 20;
md_update.Compute_Hessian_GEVP(num_evals, oversampling);
[z_update_mean, z_update_samples] = md_update.Posterior_Update_Samples();

figure;
hold on;
plot(t(2:end), z_lofi, 'color', 'black', 'LineWidth', 3);
plot(t(2:end), z_hifi, 'color', 'cyan', 'LineWidth', 3);
plot(t(2:end), z_update_mean, '--', 'color', 'red', 'LineWidth', 3);
for k = 1:num_post_samples
    plot(t(2:end), z_update_samples(:, k), 'color', [.9, .9, .9], 'LineWidth', 3);
end
plot(t(2:end), md_update.z_opt, 'color', 'black', 'LineWidth', 3);
plot(t(2:end), z_hifi, 'color', 'cyan', 'LineWidth', 3);
plot(t(2:end), z_update_mean, '--', 'color', 'red', 'LineWidth', 3);
legend({'Low-fidelity control', 'High-fidelity control', 'Update'});
set(gca, 'fontsize', 18);

%%
Jhat_update_samples = zeros(num_post_samples, 1);
for k = 1:num_post_samples
    Jhat_update_samples(k) = opt_hifi.Jhat(z_update_samples(:, k));
end
Jhat_hifi = opt_hifi.Jhat(z_hifi);
Jhat_lofi = opt_hifi.Jhat(z_lofi);
Jhat_update = opt_hifi.Jhat(z_update_mean);

save('HDSA_Results.mat');
