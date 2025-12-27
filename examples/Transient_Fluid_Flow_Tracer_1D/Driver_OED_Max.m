% Get the OED Setup ready
% OED_Setup;

% % Perform Offline OED Computations - USES data_interface
% md_oed = MD_OED_DeltaCov(opt_prob_interface, data_interface, u_prior_interface, z_prior_interface, md_hessian_analysis, oed_interface);
% md_oed.Offline_Computation();

% num_points = 12: 7.0989
% num_points = 7: 7.6511
% num_points = 5: 8.0575

num_points = 50;
Z = [z_lofi z_hifi z_hifi + 0.01 * randn(size(z_hifi, 1), num_points - 2) / norm(z_hifi)];
D = Evaluate_Discrepancy(con_hifi, con_lofi, Z);
data_interface.Set_Z_and_D(Z, D);

tdisp("Discrepancy data computed...");
% Perform Posterior Sampling (TODO: Reuse data)
md_post_sampling = MD_Posterior_Sampling(data_interface, u_prior_interface, z_prior_interface);
md_post_sampling.Compute_Posterior_Data(alpha_d, 1);

tdisp("Posterior data computed...");
% Obtain Optimal Solution Update via Continuation
num_continuation_steps = 3;
md_cont_update = MD_Continuation_Update(md_post_sampling, md_hessian_analysis, num_continuation_steps);
[~, z_cont, ~] = md_cont_update.Posterior_Update_Mean_PC_beta();
z_bar = z_cont(:, end);
tdisp("Continuation completed...");
disp(Jhat_hifi_fn(z_bar));
% 5.7337 (30 runs, alpha_u = 1)
