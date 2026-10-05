%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Import the OED setup and previously generated OED samples
OED_Setup;
oed_results = load('oed-results.mat');
num_data_points = 2;
data_interface.Set_Z_and_D(oed_results.Z_oed(:, 1:num_data_points), oed_results.D_oed(:, 1:num_data_points));

% Posterior sampling
num_samples = 10;
md_post_sampling = MD_Posterior_Sampling(data_interface, u_prior_interface, z_prior_interface);
md_post_sampling.Compute_Posterior_Data(alpha_d, num_samples, true);

% Linearization
fprintf('Linearization Update...\n');
md_update = MD_Update(md_post_sampling, md_hessian_analysis);
fprintf('Linearization Sampling...\n');
[z_lin_mean, z_lin_samples] = md_update.Posterior_Update_Samples();

% Continuation
fprintf('Continuation Update...\n');
num_continuation_steps = 3;
md_cont_update = MD_Continuation_Update(md_post_sampling, md_hessian_analysis, num_continuation_steps);
[u_mean, z_mean, beta_mean] = md_cont_update.Posterior_Update_Mean(); %#ok<NASGU>
fprintf('Continuation Sampling...\n');
tic;
[u_samples, z_samples, beta_samples] = md_cont_update.Posterior_Update_Samples(); %#ok<NASGU>
toc;

figure;
hold on;

% Calculate percentiles
plot_ptile = 100;
[z_low, z_high] = deal(prctile(z_samples, 100 - plot_ptile, 2), prctile(z_samples, plot_ptile, 2));
[z_lin_low, z_lin_high] = deal(prctile(z_lin_samples, 100 - plot_ptile, 2), prctile(z_lin_samples, plot_ptile, 2));

% Shaded regions
x_plot = 1:numel(z_mean);
plot(x_plot, z_lofi, 'b-', 'LineWidth', 2, 'DisplayName', 'LoFi');
fill([x_plot, fliplr(x_plot)], [z_lin_low(:).', fliplr(z_lin_high(:).')], ...
     [0.25 0.25 0.25], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'DisplayName', 'Linearized samples');
plot(x_plot, z_lin_mean, 'k--', 'LineWidth', 2, 'DisplayName', 'Linearized Update');
fill([x_plot, fliplr(x_plot)], [z_low(:).', fliplr(z_high(:).')], ...
     [0.25 0.75 0.25], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'DisplayName', 'Continuation samples');
plot(x_plot, z_mean, '--', 'Color', [0.25 0.75 0.25], 'LineWidth', 2, 'DisplayName', 'Continuation Update');
plot(x_plot, z_hifi, 'r-', 'LineWidth', 2, 'DisplayName', 'HiFi');

legend('Location', 'northeastoutside');
xlabel('Control index', 'Interpreter', 'latex');
ylabel('Injection/production rate', 'Interpreter', 'latex');
title('Posterior control samples for subsurface flow', 'Interpreter', 'latex');
hold off;
