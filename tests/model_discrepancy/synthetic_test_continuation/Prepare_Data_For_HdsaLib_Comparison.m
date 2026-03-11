%%
clear;
close all;
addpath(genpath('../../../src'));
rng(121234);

random_numbers = randn(10^5, 1);
writematrix(random_numbers, 'random_numbers.txt');

rng(121234);

m = 51;
x = linspace(0, 1, m)';

data_interface = MD_Data_Interface_synthetic_test_with_hessian_gevp(m);
data_interface.Load_Data();

u_prior_interface = MD_u_Prior_Interface_synthetic_test_with_hessian_gevp(m);
z_prior_interface = MD_z_Prior_Interface_synthetic_test_with_hessian_gevp(m);
md_prior_sampling = MD_Prior_Sampling(data_interface, u_prior_interface, z_prior_interface);

md_post_sampling = MD_Posterior_Sampling(data_interface, u_prior_interface, z_prior_interface);
alpha_d = 1.e-5;
num_post_samples = 1;
md_post_sampling.Compute_Posterior_Data(alpha_d, num_post_samples);

opt_prob_interface = MD_Opt_Prob_Interface_synthetic_test_with_hessian_gevp(m);
md_hessian_analysis = MD_Hessian_Analysis(opt_prob_interface, z_prior_interface);

num_evals = 20;
oversampling = 10;
md_hessian_analysis.Compute_Hessian_GEVP(data_interface.z_opt, num_evals, oversampling);

num_continuation_steps = 3;
md_cont_update = MD_Continuation_Update(md_post_sampling, md_hessian_analysis, num_continuation_steps);
[u_cont, z_cont, betas_cont] = md_cont_update.Posterior_Update_Mean_PC_beta();
z_bar = z_cont(:, end);
disp(norm(z_bar));

% save('Sabl_Output.mat', 'prior_delta', 'prior_delta_z_opt', 'post_delta_mean', 'post_delta_samples', 'post_z_mean', 'post_z_samples');
save('reference_solution.mat', 'u_cont', 'z_cont', 'betas_cont');
