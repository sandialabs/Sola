%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Import the OED
OED_Setup;
data_interface.Set_Z_and_D(z_lofi, Evaluate_Discrepancy(con_hifi, con_lofi, z_lofi))

% Posterior Sampling
num_samples = 100;
md_post_sampling = MD_Posterior_Sampling(data_interface, u_prior_interface, z_prior_interface);
md_post_sampling.Compute_Posterior_Data(alpha_d, num_samples, false);

% Continuation
num_continuation_steps = 3;
md_cont_update = MD_Continuation_Update(md_post_sampling, md_hessian_analysis, num_continuation_steps);
[u_mean, z_mean, beta_mean] = md_cont_update.Posterior_Update_Mean();
tic;
[u_samples, z_samples, beta_samples] = md_cont_update.Posterior_Update_Samples();
toc;

% Plot samples
figure;
plot(z_lofi, "b-", "LineWidth", 2, "DisplayName", "LoFi")
hold on;
plot(z_samples, "Color", [0.5 0.5 0.5], "LineWidth", 1, "HandleVisibility", "off")
plot(z_mean, "r-", "LineWidth", 2, "DisplayName", "Post-Mean")
legend("Location", "northeastoutside");
plot(z_hifi, "k-", "LineWidth", 2, "DisplayName", "HiFi")
hold off;
