%% Set up
clear;
close all;
clc;
addpath(genpath('../../src'));
load OED_Ensemble_Results.mat;

num_post_samps = 1;

M = 100;
Omega = randn(length(z_lofi), M);
z = z_lofi + 2 * z_prior_interface.Apply_W_z_Inverse_Factor(Omega);
d = con_hifi.State_Solve(z) - con_lofi.State_Solve(z);

delta_mean = cell(p, samps_per_N);
delta_samples = cell(p, samps_per_N);
for k = 1:p
    for i = 1:samps_per_N
        data_interface = MD_Data_Interface_Diff(k, i, 'OED');
        data_interface.Load_Data();
        md_update = MD_Update(opt_prob_interface, data_interface, u_prior_interface, z_prior_interface, md_hessian_analysis);
        md_update.Compute_Posterior_Data(alpha_d, num_post_samps);

        [delta_mean{k, i}, delta_samples{k, i}] = md_update.Posterior_Discrepancy_Samples(z);
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
