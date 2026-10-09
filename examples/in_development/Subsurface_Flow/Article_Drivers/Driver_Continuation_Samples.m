%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Import the OED setup and previously generated OED samples
OED_Setup;
oed_results = load('oed-results.mat');
num_data_points = 3;
data_interface.Set_Z_and_D(oed_results.Z_oed(:, 1:num_data_points), oed_results.D_oed(:, 1:num_data_points));

% Posterior sampling
num_samples = 100;
use_hybrid = false;
md_post_sampling = MD_Posterior_Sampling(data_interface, u_prior_interface, z_prior_interface);
md_post_sampling.Compute_Posterior_Data(alpha_d, num_samples, true);

% Linearization
fprintf('Linearization Update...\n');
md_update = MD_Update(md_post_sampling, md_hessian_analysis);
fprintf('Linearization Sampling...\n');
[z_lin_mean, z_lin_samples] = md_update.Posterior_Update_Samples();


% Continuation
num_continuation_steps = 7;
if ~use_hybrid
     fprintf('Continuation Update...\n');
     md_cont_update = MD_Continuation_Update(md_post_sampling, md_hessian_analysis, num_continuation_steps);
     [u_mean, z_mean, beta_mean] = md_cont_update.Posterior_Update_Mean();
     disp("Continuation Sampling...")
     % tic;
     [u_samples, z_samples, beta_samples] = md_cont_update.Posterior_Update_Samples();
     % toc;
else
     % Hybrid: continuation for the mean, linearize about it for the samples
     disp("Hybrid Update...")
     md_hybrid_update = MD_Hybrid_Continuation_Update(md_post_sampling, md_hessian_analysis, num_continuation_steps);
     [u_mean, z_mean, beta_mean] = md_hybrid_update.Posterior_Update_Mean();
     disp("Hybrid Sampling...")
     tic;
     [u_samples, z_samples, beta_samples] = md_hybrid_update.Posterior_Update_Samples();
     toc;
end

% Calculate Percentiles
plot_ptile = 90;
[z_lin_low, z_lin_high]  = deal(prctile(z_lin_samples, 50-plot_ptile/2, 2), prctile(z_lin_samples, 50+plot_ptile/2, 2));
[z_low, z_high]  = deal(prctile(z_samples, 50-plot_ptile/2, 2), prctile(z_samples, 50+plot_ptile/2, 2));

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

legend('Location', 'northeastoutside');
xlabel('Control index', 'Interpreter', 'latex');
ylabel('Injection/production rate', 'Interpreter', 'latex');
title('Posterior control samples for subsurface flow', 'Interpreter', 'latex');
ylim([-40 60])
hold off;

use_hybrid=true;
if ~use_hybrid
     num_continuation_steps = 15;
     disp("Hybrid Update...")
     md_hybrid_update = MD_Hybrid_Continuation_Update(md_post_sampling, md_hessian_analysis, num_continuation_steps);
     [u_hyb_mean, z_hyb_mean, beta_hyb_mean] = md_hybrid_update.Posterior_Update_Mean();
     disp("Hybrid Sampling...")
     tic;
     [u_hyb_samples, z_hyb_samples, beta_hyb_samples] = md_hybrid_update.Posterior_Update_Samples();
     toc;

     [z_hyb_low, z_hyb_high]  = deal(prctile(z_hyb_samples, 50-plot_ptile/2, 2), prctile(z_hyb_samples, 50+plot_ptile/2, 2));

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
     fill([x_plot, fliplr(x_plot)], [z_hyb_low(:).', fliplr(z_hyb_high(:).')], ...
          [0.75 0.25 0.75], 'FaceAlpha', 0.5, 'EdgeColor', 'none', 'DisplayName', 'Hybrid samples');
     plot(x_plot, z_mean, '--', 'Color', [0.25 0.75 0.25], 'LineWidth', 2, 'DisplayName', 'Continuation update');
     plot(x_plot, z_hifi, 'r-', 'LineWidth', 2, 'DisplayName', 'HiFi');

     legend('Location', 'northeastoutside');
     xlabel('Control index', 'Interpreter', 'latex');
     ylabel('Injection/production rate', 'Interpreter', 'latex');
     title('Posterior control samples for subsurface flow', 'Interpreter', 'latex');
     ylim([-40 60])
     hold off;
end