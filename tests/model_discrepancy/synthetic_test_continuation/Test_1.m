%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
clear;
close all;
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
num_post_samples = 10;
md_post_sampling.Compute_Posterior_Data(alpha_d, num_post_samples);

opt_prob_interface = MD_Opt_Prob_Interface_synthetic_test_continuation(m);
md_hessian_analysis = MD_Hessian_Analysis(opt_prob_interface, z_prior_interface);

md_hessian_analysis.Compute_Hessian_GEVP(data_interface.z_opt, 30, 10);

%%
num_continuation_steps = 3;
md_cont_update = MD_Continuation_Update(md_post_sampling, md_hessian_analysis, num_continuation_steps);
[u_cont, z_cont, beta_cont] = md_cont_update.Posterior_Update_Mean();
[u_ks, z_ks, beta_ks] = md_cont_update.Posterior_Update_Samples();

ref_sol = load('reference_solution.mat');
rel_err = @(val_est, val_true) norm(val_est - val_true) / norm(val_true);
u_cont_err = rel_err(u_cont, ref_sol.u_cont);
z_cont_err = rel_err(z_cont, ref_sol.z_cont);
beta_cont_err = rel_err(beta_cont, ref_sol.beta_cont);
u_ks_err = rel_err(u_ks, ref_sol.u_ks);
z_ks_err = rel_err(z_ks, ref_sol.z_ks);
beta_ks_err = rel_err(beta_ks, ref_sol.beta_ks);

ref_diff = max([u_cont_err, z_cont_err, beta_cont_err, u_ks_err, z_ks_err, beta_ks_err]);

if ref_diff > 1.e-9
    fprintf(2, '\nmodel_discrepancy/synthetic_test_continuation Test 1 failed.\n');
else
    fprintf(1, '\nmodel_discrepancy/synthetic_test_continuation Test 1 passed.\n');
end
