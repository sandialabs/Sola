%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath(genpath('../Article_B_Drivers'));

% Import the OED
OED_Setup;
oed_results = load("oed-results.mat");
num_data_points = 3;
data_interface.Set_Z_and_D(oed_results.Z_oed(:, 1:num_data_points), oed_results.D_oed(:, 1:num_data_points))

% Posterior Sampling
num_samples = 25;
md_post_sampling = MD_Posterior_Sampling(data_interface, u_prior_interface, z_prior_interface);
md_post_sampling.Compute_Posterior_Data(alpha_d, num_samples, true);

% Linearization
disp("Linearization Update...")
md_update = MD_Update(md_post_sampling, md_hessian_analysis);
disp("Linearization Sampling...")
tic;
[z_lin_mean, z_lin_samples] = md_update.Posterior_Update_Samples();
toc; 

% Continuation
disp("Continuation Update...")
num_continuation_steps = 3;
md_cont_update = MD_Continuation_Update(md_post_sampling, md_hessian_analysis, num_continuation_steps);
[u_mean, z_mean, beta_mean] = md_cont_update.Posterior_Update_Mean();
disp("Continuation Sampling...")
tic;
[u_samples, z_samples, beta_samples] = md_cont_update.Posterior_Update_Samples();
toc;

% Hybrid: continuation for the mean, linearize about it for the samples
disp("Hybrid Update...")
md_hybrid_update = MD_Hybrid_Continuation_Update(md_post_sampling, md_hessian_analysis, num_continuation_steps);
[u_hyb_mean, z_hyb_mean, beta_hyb_mean] = md_hybrid_update.Posterior_Update_Mean();
disp("Hybrid Sampling...")
tic;
[u_hyb_samples, z_hyb_samples, beta_hyb_samples] = md_hybrid_update.Posterior_Update_Samples();
toc;

% Sample-wise errors relative to continuation sampling
% rel_err = @(A) vecnorm(A - z_samples) ./ vecnorm(z_samples - z_lofi);
% fprintf("Mean rel. err. vs continuation: hybrid %.3e, linearized %.3e\n", ...
%     norm(z_hyb_mean - z_mean) / norm(z_mean - z_lofi), norm(z_lin_mean - z_mean) / norm(z_mean - z_lofi));
% fprintf("Median sample rel. err. vs continuation: hybrid %.3e, linearized %.3e\n", ...
%     median(rel_err(z_hyb_samples)), median(rel_err(z_lin_samples)));


% % Plot samples
% figure;
% plot(z_lofi, "b-", "LineWidth", 2, "DisplayName", "LoFi")
% hold on;
% plot(z_samples, "Color", [0.5 0.5 0.5], "LineWidth", 1, "HandleVisibility", "off")
% plot(z_mean, "r-", "LineWidth", 2, "DisplayName", "Post-Mean")
% plot(z_lin_samples, "Color", [0.2 0.5 0.5], "LineWidth", 1, "HandleVisibility", "off")
% plot(z_lin_mean, "m--", "LineWidth", 2, "DisplayName", "Post-Mean-Lin")
% legend("Location", "northeastoutside");
% plot(z_hifi, "k-", "LineWidth", 2, "DisplayName", "HiFi")
% hold off;

figure;
hold on;

% Calculate Percentiles
plot_ptile = 80;
[z_low, z_high]  = deal(prctile(z_samples, 100-plot_ptile, 2), prctile(z_samples, plot_ptile, 2));
[z_lin_low, z_lin_high]  = deal(prctile(z_lin_samples, 100-plot_ptile, 2), prctile(z_lin_samples, plot_ptile, 2));
[z_hyb_low, z_hyb_high] = deal(prctile(z_hyb_samples, 100-plot_ptile, 2), prctile(z_hyb_samples, plot_ptile, 2));

% Shaded regions
x = 1:numel(z_mean);
plot(x, z_lofi, "b-", "LineWidth", 2, "DisplayName", "LoFi");
fill([x, fliplr(x)],[z_lin_low(:).', fliplr(z_lin_high(:).')], ...
     [0.25 0.25 0.25], "FaceAlpha", 0.5, "EdgeColor", "none", "DisplayName", "Linearized samples");
plot(x, z_lin_mean, "k--", "LineWidth", 2, "DisplayName", "Linearized Update");
fill([x, fliplr(x)], [z_low(:).', fliplr(z_high(:).')], ...
     [0.25 0.75 0.25], "FaceAlpha", 0.5, "EdgeColor", "none", "DisplayName", "Continuation samples");
plot(x, z_mean, "--", "Color", [0.25 0.75 0.25], "LineWidth", 2, "DisplayName", "Continuation Update");
plot(x, z_hifi, "r-", "LineWidth", 2, "DisplayName", "HiFi");

legend("Location", "northeastoutside");
hold off;

% Hybrid vs continuation samples
figure;
hold on;
plot(x, z_lofi, "b-", "LineWidth", 2, "DisplayName", "LoFi");
fill([x, fliplr(x)], [z_hyb_low(:).', fliplr(z_hyb_high(:).')], ...
     [0.75 0.25 0.75], "FaceAlpha", 0.4, "EdgeColor", "none", "DisplayName", "Hybrid samples");
fill([x, fliplr(x)], [z_low(:).', fliplr(z_high(:).')], ...
     [0.25 0.75 0.25], "FaceAlpha", 0.4, "EdgeColor", "none", "DisplayName", "Continuation samples");
plot(x, z_mean, "--", "Color", [0.25 0.75 0.25], "LineWidth", 2, "DisplayName", "Continuation Update");
plot(x, z_hifi, "r-", "LineWidth", 2, "DisplayName", "HiFi");
legend("Location", "northeastoutside");
hold off;

% Very basic mean subtraction
figure;
hold on;

z_lin_low2 = z_lin_low + (z_mean - z_lin_mean);
z_lin_high2 = z_lin_high + (z_mean - z_lin_mean);

plot(x, z_lofi, "b-", "LineWidth", 2, "DisplayName", "LoFi");
fill([x, fliplr(x)],[z_lin_low2(:).', fliplr(z_lin_high2(:).')], ...
     [0.25 0.25 0.25], "FaceAlpha", 0.5, "EdgeColor", "none", "DisplayName", "Linearized samples");
plot(x, z_lin_mean, "k--", "LineWidth", 2, "DisplayName", "Linearized Update");
fill([x, fliplr(x)], [z_low(:).', fliplr(z_high(:).')], ...
     [0.25 0.75 0.25], "FaceAlpha", 0.5, "EdgeColor", "none", "DisplayName", "Continuation samples");
plot(x, z_mean, "--", "Color", [0.25 0.75 0.25], "LineWidth", 2, "DisplayName", "Continuation Update");
plot(x, z_hifi, "r-", "LineWidth", 2, "DisplayName", "HiFi");

legend("Location", "northeastoutside");
hold off;


fprintf("Mean rel. err. vs continuation: hybrid %.3e, linearized %.3e\n", ...
    norm(z_hyb_mean - z_mean) / norm(z_mean - z_lofi), norm(z_lin_mean - z_mean) / norm(z_mean - z_lofi));
fprintf("Low-end sample rel. err. vs continuation: hybrid %.3e, linearized %.3e\n", ...
    norm(z_hyb_low - z_low)/norm(z_low), norm(z_lin_low - z_low)/norm(z_low));
fprintf("High-end sample rel. err. vs continuation: hybrid %.3e, linearized %.3e\n", ...
    norm(z_hyb_high - z_high)/norm(z_high), norm(z_lin_low - z_low)/norm(z_low));
