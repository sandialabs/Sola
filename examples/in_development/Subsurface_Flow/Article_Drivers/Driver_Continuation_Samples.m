%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Import the OED setup and previously generated OED samples
OED_Setup;
oed_results = load('oed-results.mat');
num_data_points = 5;
data_interface.Set_Z_and_D(oed_results.Z_oed(:, 1:num_data_points), oed_results.D_oed(:, 1:num_data_points));

% Posterior sampling
num_samples = 50;
md_post_sampling = MD_Posterior_Sampling(data_interface, u_prior_interface, z_prior_interface);
md_post_sampling.Compute_Posterior_Data(alpha_d, num_samples, true);

% Linearization
fprintf('Linearization Update...\n');
md_update = MD_Update(md_post_sampling, md_hessian_analysis);
fprintf('Linearization Sampling...\n');
[z_lin_mean, z_lin_samples] = md_update.Posterior_Update_Samples();

% Continuation
fprintf('Continuation Update...\n');
md_cont_update = MD_Continuation_Update(md_post_sampling, md_hessian_analysis);

mean_options = struct();
mean_options.num_steps = 7;
% mean_options.time_step.use = true;
[u_mean, z_mean, beta_mean] = md_cont_update.Posterior_Update_Mean(mean_options);


sample_options = struct();
sample_options.num_mean_to_sample_steps = 4;
% sample_options.time_step.use = true;
sample_options.materialize_breve_samples = true;
sample_options.use_parallel_sampling = true;

pool = parpool("Threads");
disp("Continuation Sampling...")
tic;
[u_samples, z_samples, beta_samples] = md_cont_update.Posterior_Update_Samples_From_Mean(u_mean, z_mean, beta_mean, sample_options);
toc;
delete(pool)

% Calculate Percentiles
plot_ptile = 90;
get_percentile = @(samples, ptile) deal(prctile(samples, 50-ptile/2, 2), prctile(samples, 50+ptile/2, 2));
[z_lin_low, z_lin_high] = get_percentile(z_lin_samples, plot_ptile);
[z_low, z_high] = get_percentile(z_samples, plot_ptile);

figure;
hold on;

% Shaded regions
x_plot = 1:numel(z_mean);
plot(x_plot, z_lofi, 'b-', 'LineWidth', 2, 'DisplayName', 'LoFi');
fill([x_plot, fliplr(x_plot)], [z_lin_low(:).', fliplr(z_lin_high(:).')], ...
     [0.25 0.25 0.25], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'DisplayName', 'Linearized samples');
plot(x_plot, z_lin_mean, 'k--', 'LineWidth', 2, 'DisplayName', 'Linearized update');
fill([x_plot, fliplr(x_plot)], [z_low(:).', fliplr(z_high(:).')], ...
     [0.25 0.75 0.25], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'DisplayName', 'Continuation samples');
plot(x_plot, z_mean, '--', 'Color', [0.25 0.75 0.25], 'LineWidth', 2, 'DisplayName', 'Continuation update');
plot(x_plot, z_hifi, 'r-', 'LineWidth', 2, 'DisplayName', 'HiFi');
plot(x_plot, z_best_proj, 'c-', 'LineWidth', 2, 'DisplayName', 'HiFi Proj');

legend('Location', 'northeastoutside');
xlabel('Control index', 'Interpreter', 'latex');
ylabel('Injection/production rate', 'Interpreter', 'latex');
title('Posterior control samples for subsurface flow', 'Interpreter', 'latex');
ylim([-40 60])
hold off;

% exportgraphics(gca, ['new_plot_', num2str(num_data_points), '.pdf'])