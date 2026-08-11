%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
clear;
close all;
rng(121234);

random_numbers = randn(10^5, 1);
writematrix(random_numbers, 'random_numbers.txt');

rng(121234);

m = 51;
x = linspace(0, 1, m)';

data_interface = MD_Data_Interface_synthetic_test_continuation(m);
data_interface.Load_Data();

u_prior_interface = MD_u_Prior_Interface_synthetic_test_continuation(m);
z_prior_interface = MD_z_Prior_Interface_synthetic_test_continuation(m);
md_prior_sampling = MD_Prior_Sampling(data_interface, u_prior_interface, z_prior_interface);

md_post_sampling = MD_Posterior_Sampling(data_interface, u_prior_interface, z_prior_interface);
alpha_d = 1.e-5;
num_post_samples = 100;
md_post_sampling.Compute_Posterior_Data(alpha_d, num_post_samples);

opt_prob_interface = MD_Opt_Prob_Interface_synthetic_test_continuation(m);
md_hessian_analysis = MD_Hessian_Analysis(opt_prob_interface, z_prior_interface);
num_evals = 20;
oversampling = 10;
md_hessian_analysis.Compute_Hessian_GEVP(data_interface.z_opt, num_evals, oversampling);

%%
num_continuation_steps = 3;
md_cont_update = MD_Continuation_Update(md_post_sampling, md_hessian_analysis, num_continuation_steps);
[u_cont, z_cont, beta_cont] = md_cont_update.Posterior_Update_Mean();
% [u_ks, z_ks, beta_ks] = md_cont_update.Posterior_Update_Samples();
u_k = u_cont(:, end);
z_k = z_cont(:, end);
beta_k = beta_cont(:, end);
disp(norm(z_cont(:, end)));

save('Sabl_output.mat', 'u_k', 'z_k', 'beta_k');
save('reference_solution.mat', 'u_cont', 'z_cont', 'beta_cont');
%  'u_ks', 'z_ks', 'beta_ks'

% ------------------------------------------------------------
% Printers
% ------------------------------------------------------------

low_fidelity_state = @(z_in) z_in.^3;
high_fidelity_state = @(z_in) 1.2 * z_in.^3;
z_lf_opt = 1 + x;
z_hf_opt = 1.2^(-1/3) * z_lf_opt;

u_lf_at_lf_opt  = low_fidelity_state(z_lf_opt);
u_lf_at_updated = low_fidelity_state(z_k);
u_hf_at_lf_opt  = high_fidelity_state(z_lf_opt);
u_hf_at_updated = high_fidelity_state(z_k);
u_hf_at_hf_opt = high_fidelity_state(z_hf_opt);

J_lf_at_lf_opt  = opt_prob_interface.Objective_Function(u_lf_at_lf_opt,  z_lf_opt);
J_lf_at_updated = opt_prob_interface.Objective_Function(u_lf_at_updated, z_k);
J_hf_at_lf_opt  = opt_prob_interface.Objective_Function(u_hf_at_lf_opt,  z_lf_opt);
J_hf_at_updated = opt_prob_interface.Objective_Function(u_hf_at_updated, z_k);
J_hf_at_hf_opt = opt_prob_interface.Objective_Function(u_hf_at_hf_opt, z_hf_opt);

% ------------------------------------------------------------
% Print results
% ------------------------------------------------------------

fprintf('-----------------------------------------------------\n');
fprintf('Actual objective comparison\n');
fprintf('\nJ_LF at low-fidelity optimum:       %g\n', J_lf_at_lf_opt);
fprintf('J_LF at updated solution:           %g\n', J_lf_at_updated);
fprintf('\nJ_HF at low-fidelity optimum:       %g\n', J_hf_at_lf_opt);
fprintf('J_HF at updated solution:           %g\n', J_hf_at_updated);
fprintf('J_HF at exact HF optimum:           %g\n', J_hf_at_hf_opt);
fprintf('\nHF improvement factor:              %g%%\n', 100.0 * (1.0 - J_hf_at_updated / J_hf_at_lf_opt));
fprintf('\n||z_LF_opt - z_HF_opt||:            %g\n', norm(z_lf_opt - z_hf_opt));
fprintf('||z_updated - z_HF_opt||:           %g\n', norm(z_cont(:, end) - z_hf_opt));
