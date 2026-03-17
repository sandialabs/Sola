%% Set up
clear;
close all;
clc;
addpath(genpath('../../src'));
load Optimization_Results.mat;
rng(2451423);

T = 5;
N = 100;
obj_hifi = Mass_Spring_Objective_HiFi(T, N);
obj_lofi = Mass_Spring_Objective_LoFi(obj_hifi);
con_hifi = Mass_Spring_Coupled(T, N);
con_lofi = Mass_Spring_LoFi(con_hifi);
opt_hifi = Reduced_Space_Optimization(obj_hifi, con_hifi);
opt_lofi = Reduced_Space_Optimization(obj_lofi, con_lofi);
t = con_hifi.t_mesh;

%%
data_interface = MD_Data_Interface_Mass_Spring();
data_interface.Load_Data();

alpha_u = 50;
alpha_z = 1.e-2;
u_prior_interface = MD_Elliptic_u_Prior_Interface_Mass_Spring(alpha_u, opt_lofi);
z_prior_interface = MD_Elliptic_z_Prior_Interface_Mass_Spring(alpha_z, opt_lofi);

%%
num_prior_samples = 500;
md_prior_sampling = MD_Prior_Sampling(data_interface, u_prior_interface, z_prior_interface);

prior_delta_samples = md_prior_sampling.Prior_Discrepancy_Samples_at_z_opt(num_prior_samples);

figure;
hold on;
plot(t, prior_delta_samples(1:4:end, :), 'LineWidth', 3, 'color', [.9, .9, .9]);
plot(t, prior_delta_samples(1:4:end, 1:5), 'LineWidth', 3);
xlabel('Time');
ylabel('$x_1$', 'Interpreter', 'latex');
set(gca, 'fontsize', 18);

figure;
hold on;
plot(t, prior_delta_samples(2:4:end, :), 'LineWidth', 3, 'color', [.9, .9, .9]);
plot(t, prior_delta_samples(2:4:end, 1:5), 'LineWidth', 3);
xlabel('Time');
ylabel('$v_1$', 'Interpreter', 'latex');
set(gca, 'fontsize', 18);

figure;
hold on;
plot(t, prior_delta_samples(3:4:end, :), 'LineWidth', 3, 'color', [.9, .9, .9]);
plot(t, prior_delta_samples(3:4:end, 1:5), 'LineWidth', 3);
xlabel('Time');
ylabel('$x_2$', 'Interpreter', 'latex');
set(gca, 'fontsize', 18);

figure;
hold on;
plot(t, prior_delta_samples(4:4:end, :), 'LineWidth', 3, 'color', [.9, .9, .9]);
plot(t, prior_delta_samples(4:4:end, 1:5), 'LineWidth', 3);
xlabel('Time');
ylabel('$v_2$', 'Interpreter', 'latex');
set(gca, 'fontsize', 18);

z_samples = z_prior_interface.Sample_with_Covariance_W_z_Inverse(5);
figure;
hold on;
plot(t(2:end), z_samples, 'LineWidth', 3);
set(gca, 'fontsize', 18);

%%
md_post_sampling = MD_Posterior_Sampling(data_interface, u_prior_interface, z_prior_interface);
alpha_d = 5.e-4;
num_post_samples = 500;
md_post_sampling.Compute_Posterior_Data(alpha_d, num_post_samples);
Z_test = zeros(N - 1, 2);
Z_test(:, 1) = Z;
Z_test(:, 2) = 2 * ones(N - 1, 1);
[delta_mean, delta_samples] = md_post_sampling.Posterior_Discrepancy_Samples(Z_test);

figure;
hold on;
plot(t(2:end), Z_test(:, 1), 'LineWidth', 3);
plot(t(2:end), Z_test(:, 2), 'LineWidth', 3);
legend({'z_1', 'z_2'});
set(gca, 'fontsize', 18);

figure;
hold on;
plot(t, md_post_sampling.post_data.D(1:4:end, 1), 'color', 'black', 'LineWidth', 3);
plot(t, delta_mean{1}(1:4:end), '--', 'color', 'red', 'LineWidth', 3);
for k = 1:num_post_samples
    plot(t, delta_samples{1}(1:4:end, k), 'color', [.9, .9, .9], 'LineWidth', 3);
end
plot(t, md_post_sampling.post_data.D(1:4:end, 1), 'color', 'black', 'LineWidth', 3);
plot(t, delta_mean{1}(1:4:end), '--', 'color', 'red', 'LineWidth', 3);
xlabel('Time');
ylabel('$x_1$', 'Interpreter', 'latex');
set(gca, 'fontsize', 18);

figure;
hold on;
plot(t, md_post_sampling.post_data.D(2:4:end, 1), 'color', 'black', 'LineWidth', 3);
plot(t, delta_mean{1}(2:4:end), '--', 'color', 'red', 'LineWidth', 3);
for k = 1:num_post_samples
    plot(t, delta_samples{1}(2:4:end, k), 'color', [.9, .9, .9], 'LineWidth', 3);
end
plot(t, md_post_sampling.post_data.D(2:4:end, 1), 'color', 'black', 'LineWidth', 3);
plot(t, delta_mean{1}(2:4:end), '--', 'color', 'red', 'LineWidth', 3);
xlabel('Time');
ylabel('$v_1$', 'Interpreter', 'latex');
set(gca, 'fontsize', 18);

figure;
hold on;
plot(t, md_post_sampling.post_data.D(3:4:end, 1), 'color', 'black', 'LineWidth', 3);
plot(t, delta_mean{1}(3:4:end), '--', 'color', 'red', 'LineWidth', 3);
for k = 1:num_post_samples
    plot(t, delta_samples{1}(3:4:end, k), 'color', [.9, .9, .9], 'LineWidth', 3);
end
plot(t, md_post_sampling.post_data.D(3:4:end, 1), 'color', 'black', 'LineWidth', 3);
plot(t, delta_mean{1}(3:4:end), '--', 'color', 'red', 'LineWidth', 3);
xlabel('Time');
ylabel('$x_2$', 'Interpreter', 'latex');
set(gca, 'fontsize', 18);

figure;
hold on;
plot(t, md_post_sampling.post_data.D(4:4:end, 1), 'color', 'black', 'LineWidth', 3);
plot(t, delta_mean{1}(4:4:end), '--', 'color', 'red', 'LineWidth', 3);
for k = 1:num_post_samples
    plot(t, delta_samples{1}(4:4:end, k), 'color', [.9, .9, .9], 'LineWidth', 3);
end
plot(t, md_post_sampling.post_data.D(4:4:end, 1), 'color', 'black', 'LineWidth', 3);
plot(t, delta_mean{1}(4:4:end), '--', 'color', 'red', 'LineWidth', 3);
xlabel('Time');
ylabel('$v_2$', 'Interpreter', 'latex');
set(gca, 'fontsize', 18);

figure;
hold on;
plot(t, delta_mean{2}(1:4:end), '--', 'color', 'red', 'LineWidth', 3);
for k = 1:num_post_samples
    plot(t, delta_samples{2}(1:4:end, k), 'color', [.9, .9, .9], 'LineWidth', 3);
end
plot(t, delta_mean{2}(1:4:end), '--', 'color', 'red', 'LineWidth', 3);
xlabel('Time');
ylabel('$x_1$', 'Interpreter', 'latex');
set(gca, 'fontsize', 18);

figure;
hold on;
plot(t, delta_mean{2}(2:4:end), '--', 'color', 'red', 'LineWidth', 3);
for k = 1:num_post_samples
    plot(t, delta_samples{2}(2:4:end, k), 'color', [.9, .9, .9], 'LineWidth', 3);
end
plot(t, delta_mean{2}(2:4:end), '--', 'color', 'red', 'LineWidth', 3);
xlabel('Time');
ylabel('$v_1$', 'Interpreter', 'latex');
set(gca, 'fontsize', 18);

figure;
hold on;
plot(t, delta_mean{2}(3:4:end), '--', 'color', 'red', 'LineWidth', 3);
for k = 1:num_post_samples
    plot(t, delta_samples{2}(3:4:end, k), 'color', [.9, .9, .9], 'LineWidth', 3);
end
plot(t, delta_mean{2}(3:4:end), '--', 'color', 'red', 'LineWidth', 3);
xlabel('Time');
ylabel('$x_2$', 'Interpreter', 'latex');
set(gca, 'fontsize', 18);

figure;
hold on;
plot(t, delta_mean{2}(4:4:end), '--', 'color', 'red', 'LineWidth', 3);
for k = 1:num_post_samples
    plot(t, delta_samples{2}(4:4:end, k), 'color', [.9, .9, .9], 'LineWidth', 3);
end
plot(t, delta_mean{2}(4:4:end), '--', 'color', 'red', 'LineWidth', 3);
xlabel('Time');
ylabel('$v_2$', 'Interpreter', 'latex');
set(gca, 'fontsize', 18);

%%
opt_prob_interface = MD_Opt_Prob_Interface_Sabl(opt_lofi, data_interface);
md_hessian_analysis = MD_Hessian_Analysis(opt_prob_interface, z_prior_interface);

num_evals = 6;
oversampling = 4;
md_hessian_analysis.Compute_Hessian_GEVP(data_interface.z_opt, num_evals, oversampling);

md_update = MD_Update(md_post_sampling, md_hessian_analysis);
[z_update_mean, z_update_samples] = md_update.Posterior_Update_Samples();

figure;
hold on;
plot(t(2:end), z_tilde, 'color', 'black', 'LineWidth', 3);
plot(t(2:end), z_star, 'color', 'cyan', 'LineWidth', 3);
plot(t(2:end), z_update_mean, '--', 'color', 'red', 'LineWidth', 3);
for k = 1:num_post_samples
    plot(t(2:end), z_update_samples(:, k), 'color', [.9, .9, .9], 'LineWidth', 3);
end
plot(t(2:end), z_tilde, 'color', 'black', 'LineWidth', 3);
plot(t(2:end), z_star, 'color', 'cyan', 'LineWidth', 3);
plot(t(2:end), z_update_mean, '--', 'color', 'red', 'LineWidth', 3);
legend({'Low-fidelity control', 'High-fidelity control', 'Update'});
set(gca, 'fontsize', 18);

%%
Jhat_update_samples = zeros(num_post_samples, 1);
for k = 1:num_post_samples
    Jhat_update_samples(k) = opt_hifi.Jhat(z_update_samples(:, k));
end
Jhat_hifi = opt_hifi.Jhat(z_star);
Jhat_lofi = opt_hifi.Jhat(z_tilde);
Jhat_update = opt_hifi.Jhat(z_update_mean);

figure;
hold on;
plot([Jhat_lofi, Jhat_lofi], [0, 100], 'color', 'black', 'LineWidth', 3);
plot([Jhat_hifi, Jhat_hifi], [0, 100], 'color', 'cyan', 'LineWidth', 3);
plot([Jhat_update, Jhat_update], [0, 100], 'color', 'red', 'LineWidth', 3);
histogram(Jhat_update_samples);
plot([Jhat_lofi, Jhat_lofi], [0, 100], 'color', 'black', 'LineWidth', 3);
plot([Jhat_hifi, Jhat_hifi], [0, 100], 'color', 'cyan', 'LineWidth', 3);
plot([Jhat_update, Jhat_update], [0, 100], 'color', 'red', 'LineWidth', 3);
legend({'Low-fidelity objective', 'High-fidelity objective', 'Update objective'});
set(gca, 'fontsize', 18);

figure;
hold on;
plot([Jhat_lofi, Jhat_lofi], [0, 100], 'color', 'black', 'LineWidth', 3);
plot([Jhat_hifi, Jhat_hifi], [0, 100], 'color', 'cyan', 'LineWidth', 3);
plot([Jhat_update, Jhat_update], [0, 100], 'color', 'red', 'LineWidth', 3);
histogram(Jhat_update_samples);
plot([Jhat_lofi, Jhat_lofi], [0, 100], 'color', 'black', 'LineWidth', 3);
plot([Jhat_hifi, Jhat_hifi], [0, 100], 'color', 'cyan', 'LineWidth', 3);
plot([Jhat_update, Jhat_update], [0, 100], 'color', 'red', 'LineWidth', 3);
legend({'Low-fidelity objective', 'High-fidelity objective', 'Update objective'});
xlim([0, 0.05]);
set(gca, 'fontsize', 18);

save('HDSA_Results.mat');
