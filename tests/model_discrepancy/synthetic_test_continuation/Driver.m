%%
clear;
close all;
addpath(genpath('../../../src'));
rng(121234);

suppress_figures = true;

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

%%
num_continuation_steps = 3;
md_cont_update = MD_Continuation_Update(md_post_sampling, md_hessian_analysis, num_continuation_steps);
[u_cont, z_cont, betas_cont] = md_cont_update.Posterior_Update_Mean();
z_bar = z_cont(:, end);
beta_bar = betas_cont(:, end);

u_cont_ref = load('reference_solution.mat').u_cont;
z_cont_ref = load('reference_solution.mat').z_cont;
betas_cont_ref = load('reference_solution.mat').betas_cont;
ref_diff = max([norm(u_cont_ref - u_cont) / norm(u_cont), norm(z_cont_ref - z_cont) / norm(z_cont), norm(betas_cont_ref - betas_cont) / norm(betas_cont)]);
if ref_diff > 1.e-9
    disp('model_discrepancy_continuation difference:');
    disp(ref_diff);
end

% norm(md_cont_update.Apply_Discrepancy_z_theta_Hessian(u_cont_ref(:, end)))
% norm(md_cont_update.Apply_Discrepancy_z_Jacobian_transpose(u_cont_ref(:, end), 1.0))
% u_n = u_cont_ref(:, end);
% norm(md_cont_update.Apply_Parameterized_RS_Hessian_Inverse_beta(beta_bar, u_n, beta_bar, 0.5))
% norm(md_cont_update.Apply_Parameterized_RS_Hessian_Inverse_beta_noCG(beta_bar, u_n, beta_bar, 0.5))

% u_n = u_cont_ref(:, end);
% beta_in = beta_bar;
% beta_n = beta_bar;
% t_n = 0.5;

% I = eye(length(beta_n));
% H = I;
% for i = 1:size(beta_n, 1)
%     H(:, i) = md_cont_update.Apply_Parameterized_RS_Hessian_beta(I(:, i), u_n, beta_n, t_n);
% end
% beta_out_full = H \ beta_in;

% [beta_out, flag, relres, iter, resvec] = pcg(@(x) H * x, beta_in, 1e-7, length(beta_n)+10);
% if flag ~= 0
%     disp('CG did not converge');
%     disp(flag);
% end
% disp(iter);
% disp(relres);
% disp(100*norm(beta_out - beta_out_full)/norm(beta_out_full))
