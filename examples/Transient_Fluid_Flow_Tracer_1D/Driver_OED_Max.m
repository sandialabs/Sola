Z = [z_lofi z_hifi z_hifi + 0.01 * randn(size(z_hifi, 1), 30) / norm(z_hifi)];
D = Evaluate_Discrepancy(con_hifi, con_lofi, Z);
data_interface.Set_Z_and_D(Z, D);

disp("Discrepancy data computed...");
% Perform Posterior Sampling (TODO: Reuse data)
md_post_sampling = MD_Posterior_Sampling(data_interface, u_prior_interface, z_prior_interface);
md_post_sampling.Compute_Posterior_Data(alpha_d, 1);

disp("Posterior data computed...");
% Obtain Optimal Solution Update via Continuation
num_continuation_steps = 3;
md_cont_update = MD_Continuation_Update(md_post_sampling, md_hessian_analysis, num_continuation_steps);
[~, z_cont, ~] = md_cont_update.Posterior_Update_Mean_PC_beta();
z_bar = z_cont(:, end);
disp(Jhat_hifi_fn(z_bar));
% 5.7337
