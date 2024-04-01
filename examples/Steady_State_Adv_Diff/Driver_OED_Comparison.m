%% Set up
clear;
close all;
clc;
addpath(genpath('../../src'));
load OED_Ensemble_Results.mat;

num_post_samps = 1;

Jhat_lofi = opt_hifi.Jhat(z_lofi);
Jhat_hifi = opt_hifi.Jhat(z_hifi);

z_error_0 = sqrt((z_lofi - z_hifi)' * z_prior_interface.Apply_M_z(z_lofi - z_hifi)) / sqrt(z_hifi' * z_prior_interface.Apply_M_z(z_hifi));

md_post_sampling = MD_Posterior_Sampling(data_interface, u_prior_interface, z_prior_interface);
md_post_sampling.Compute_Posterior_Data(alpha_d, num_post_samps);

md_update = MD_Update(md_post_sampling, md_hessian_analysis);
[z_update_mean, z_update_samples] = md_update.Posterior_Update_Samples();
Jhat_oed_1 = opt_hifi.Jhat(z_update_mean);
z_error_1 = sqrt((z_update_mean - z_hifi)' * z_prior_interface.Apply_M_z(z_update_mean - z_hifi)) / sqrt(z_hifi' * z_prior_interface.Apply_M_z(z_hifi));

Jhat_oed = zeros(p, samps_per_N);
oed_z_error = zeros(p, samps_per_N);
Jhat_rand = zeros(p, samps_per_N);
rand_z_error = zeros(p, samps_per_N);
for k = 1:p
    for i = 1:samps_per_N
        data_interface = MD_Data_Interface_Diff(k, i, 'OED');
        data_interface.Load_Data();
        md_post_sampling = MD_Posterior_Sampling(data_interface, u_prior_interface, z_prior_interface);
        md_post_sampling.Compute_Posterior_Data(alpha_d, num_post_samps);
        md_update = MD_Update(md_post_sampling, md_hessian_analysis);
        [z_update_mean, z_update_samples] = md_update.Posterior_Update_Samples();
        Jhat_oed(k, i) = opt_hifi.Jhat(z_update_mean);
        oed_z_error(k, i) = sqrt((z_update_mean - z_hifi)' * z_prior_interface.Apply_M_z(z_update_mean - z_hifi)) / sqrt(z_hifi' * z_prior_interface.Apply_M_z(z_hifi));

        data_interface = MD_Data_Interface_Diff(k, i, 'Random');
        data_interface.Load_Data();
        md_post_sampling = MD_Posterior_Sampling(data_interface, u_prior_interface, z_prior_interface);
        md_post_sampling.Compute_Posterior_Data(alpha_d, num_post_samps);
        md_update = MD_Update(md_post_sampling, md_hessian_analysis);
        [z_update_mean, z_update_samples] = md_update.Posterior_Update_Samples();
        Jhat_rand(k, i) = opt_hifi.Jhat(z_update_mean);
        rand_z_error(k, i) = sqrt((z_update_mean - z_hifi)' * z_prior_interface.Apply_M_z(z_update_mean - z_hifi)) / sqrt(z_hifi' * z_prior_interface.Apply_M_z(z_hifi));
    end
end

colors = zeros(7, 3);
colors(1, :) = [0 0.4470 0.7410];
colors(2, :) = [0.8500 0.3250 0.0980];
colors(3, :) = [0.9290 0.6940 0.1250];
colors(4, :) = [0.4940 0.1840 0.5560];
colors(5, :) = [0.4660 0.6740 0.1880];
colors(6, :) = [0.3010 0.7450 0.9330];
colors(7, :) = [0.6350 0.0780 0.1840];

mark_size = 10;
figure;
hold on;
plot(0, Jhat_lofi, '.', 'MarkerSize', 3 * mark_size, 'color', colors(1, :));
plot(0, Jhat_hifi, '.', 'MarkerSize', 3 * mark_size, 'color', colors(2, :));
plot(1, Jhat_oed_1, '.', 'MarkerSize', 3 * mark_size, 'color', colors(3, :));
legend({'Low-fidelity', 'High-fidelity', 'One evaluation update'}, 'FontSize', 18);
for k = 1:p
    plot(N_range(k) - 0.1, Jhat_oed(k, :), 'o', 'MarkerSize', mark_size, 'color', colors(4, :), 'HandleVisibility', 'off');
    plot(N_range(k) - 0.1, mean(Jhat_oed(k, :)), 'd', 'MarkerSize', mark_size, 'color', colors(5, :), 'HandleVisibility', 'off');
    plot(N_range(k) + 0.1, Jhat_rand(k, :), 'x', 'MarkerSize', mark_size, 'color', colors(6, :), 'HandleVisibility', 'off');
    plot(N_range(k) + 0.1, mean(Jhat_rand(k, :)), 'd', 'MarkerSize', mark_size, 'color', colors(7, :), 'HandleVisibility', 'off');
end

mark_size = 10;
figure;
hold on;
plot(0, z_error_0, '.', 'MarkerSize', 3 * mark_size, 'color', colors(1, :));
plot(1, z_error_1, '.', 'MarkerSize', 3 * mark_size, 'color', colors(3, :));
for k = 1:p
    plot(N_range(k) - 0.1, oed_z_error(k, :), 'o', 'MarkerSize', mark_size, 'color', colors(4, :), 'HandleVisibility', 'off');
    plot(N_range(k) - 0.1, mean(oed_z_error(k, :)), 'd', 'MarkerSize', mark_size, 'color', colors(5, :), 'HandleVisibility', 'off');
    plot(N_range(k) + 0.1, rand_z_error(k, :), 'x', 'MarkerSize', mark_size, 'color', colors(6, :), 'HandleVisibility', 'off');
    plot(N_range(k) + 0.1, mean(rand_z_error(k, :)), 'd', 'MarkerSize', mark_size, 'color', colors(7, :), 'HandleVisibility', 'off');
end
