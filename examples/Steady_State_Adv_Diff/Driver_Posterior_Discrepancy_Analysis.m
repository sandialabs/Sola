%% Set up
clear;
close all;
clc;
addpath(genpath('../../src'));
load OED_Ensemble_Results.mat;

num_post_samps = 100;

M = 100;
z = z_lofi + 40 * z_prior_interface.Sample_with_Covariance_W_z_Inverse(M);
d = con_hifi.State_Solve(z) - con_lofi.State_Solve(z);

delta_mean = cell(p, samps_per_N);
delta_samples = cell(p, samps_per_N);
for k = 1:p
    for i = 1:samps_per_N
        data_interface = MD_Data_Interface_Diff(k, i, 'OED');
        data_interface.Load_Data();
        md_post_sampling = MD_Posterior_Sampling(data_interface, u_prior_interface, z_prior_interface);
        md_post_sampling.Compute_Posterior_Data(alpha_d, num_post_samps);

        [delta_mean{k, i}, delta_samples{k, i}] = md_post_sampling.Posterior_Discrepancy_Samples(z);
    end
end

normalize = zeros(size(z, 2), 1);
for j = 1:size(z, 2)
    normalize(j) = sqrt((d(:, j))' * u_prior_interface.Apply_M_u(d(:, j)));
end

mean_error = zeros(p, samps_per_N, size(z, 2));
for k = 1:p
    for i = 1:samps_per_N
        for j = 1:size(z, 2)
            mean_error(k, i, j) = sqrt((delta_mean{k, i}{j} - d(:, j))' * u_prior_interface.Apply_M_u(delta_mean{k, i}{j} - d(:, j))) / normalize(j);
        end
    end
end

figure;
hold on;
for k = 1:p
    plot(N_range(k), mean(mean_error(k, :, :), 3), 'o', 'MarkerSize', 10);
end

d_hifi = con_hifi.State_Solve(z_hifi) - con_lofi.State_Solve(z_hifi);
[delta_mean_hifi, delta_samples_hifi] = md_post_sampling.Posterior_Discrepancy_Samples(z_hifi);
figure;
hold on;
plot(x, d_hifi, 'LineWidth', 3, 'Color', 'cyan');
plot(x, delta_mean_hifi{1}, 'LineWidth', 3, 'Color', 'red');
for k = 1:num_post_samps
    plot(x, delta_samples_hifi{1}(:, k), 'LineWidth', 3, 'Color', [.9, .9, .9]);
end
plot(x, d_hifi, 'LineWidth', 3, 'Color', 'cyan');
plot(x, delta_mean_hifi{1}, 'LineWidth', 3, 'Color', 'red');
ylim([-20, 10]);
legend({'Discrepancy', 'Discrepancy Approximation'});
title('Discrepancy error at high-fidelity solution');

[delta_mean_data, delta_samples_data] = md_post_sampling.Posterior_Discrepancy_Samples(md_post_sampling.data_interface.Z);
for j = 1:size(md_post_sampling.data_interface.Z, 2)
    figure;
    hold on;
    plot(x, md_post_sampling.data_interface.D(:, j), 'LineWidth', 3, 'Color', 'cyan');
    plot(x, delta_mean_data{j}, 'LineWidth', 3, 'Color', 'red');
    for k = 1:num_post_samps
        plot(x, delta_samples_data{j}(:, k), 'LineWidth', 3, 'Color', [.9, .9, .9]);
    end
    plot(x, md_post_sampling.data_interface.D(:, j), 'LineWidth', 3, 'Color', 'cyan');
    plot(x, delta_mean_data{j}, 'LineWidth', 3, 'Color', 'red');
    ylim([-20, 10]);
    legend({'Discrepancy', 'Discrepancy Approximation'});
    title(['Discrepancy error at training data ', num2str(j)]);
end

figure;
hold on;
plot(x, z_hifi, 'LineWidth', 3, 'Color', 'cyan');
plot(x, md_post_sampling.data_interface.Z(:, 1), 'LineWidth', 3, 'Color', 'red');
plot(x, md_post_sampling.data_interface.Z(:, 2), 'LineWidth', 3, 'Color', 'black');
plot(x, md_post_sampling.data_interface.Z(:, 3), 'LineWidth', 3, 'Color', 'black');
plot(x, md_post_sampling.data_interface.Z(:, 4), 'LineWidth', 3, 'Color', 'black');
legend({'High-fidelity control', 'Low-fidelity control', 'Training Data'});
