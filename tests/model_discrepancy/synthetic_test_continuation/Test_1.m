%%
clear;
close all;
rng(121234);

suppress_figures = true;

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

%%
num_continuation_steps = 3;
md_cont_update = MD_Continuation_Update(md_post_sampling, md_hessian_analysis, num_continuation_steps);
[u_cont, z_cont, beta_cont] = md_cont_update.Posterior_Update_Mean();
[u_ks, z_ks, beta_ks] = md_cont_update.Posterior_Update_Samples();
ref_sol = load('reference_solution.mat');

% disp(md_cont_update.Jhat_Posterior(z_cont, 0));
% TODO: Fix it so that those functions are exposed
rel_err = @(val_est, val_true) norm(val_est - val_true) / norm(val_true);
u_cont_err = rel_err(u_cont, ref_sol.u_cont);
z_cont_err = rel_err(z_cont, ref_sol.z_cont);
beta_cont_err = rel_err(beta_cont, ref_sol.beta_cont);
u_ks_err = rel_err(u_ks, ref_sol.u_ks);
z_ks_err = rel_err(z_ks, ref_sol.z_ks);
beta_ks_err = rel_err(beta_ks, ref_sol.beta_ks);

ref_diff = max([u_cont_err, z_cont_err, beta_cont_err, u_ks_err, z_ks_err, beta_ks_err]);

if ref_diff > 1.e-9
    fprintf(2, '\nmodel_discrepancy/synthetic_test_continuation failed.\n');
else
    fprintf(1, '\nmodel_discrepancy/synthetic_test_continuation passed.\n');
end

% ------
if false
    % Setup
    rng(1);
    z_n = data_interface.z_opt + 1;
    t_n = 1.0;
    sample_idx = 1;
    v = randn(size(data_interface.z_opt));
    w = randn(size(data_interface.u_opt));

    % Check individual sample code
    [delta_mean, delta_samples] = md_post_sampling.Posterior_Discrepancy_Samples(data_interface.z_opt);
    delta_mean = delta_mean{1};
    delta_samples = delta_samples{1};
    % disp(norm(u_out - u_out_fd)/norm(u_out))

    % Test Discrepancy evaluations
    delta_mean_1 = md_cont_update.Discrepancy_Evaluation_Mean(data_interface.z_opt, 1.0);
    delta_samples_1 = md_cont_update.Discrepancy_Evaluation_Sample(data_interface.z_opt, 1.0, 1);
    % fprintf('Rel. err in discrep @ mean:   %.3e\n', norm(delta_mean_1 - delta_mean)/norm(delta_mean))
    % fprintf('Rel. err in discrep @ sample:   %.3e\n', norm(delta_samples_1 - delta_samples(:, 1))/norm(delta_samples(:, 1)))

    % Finite Difference check for Apply_Discrepancy_z_Jacobian_Sample
    Jv = md_cont_update.Apply_Discrepancy_z_Jacobian_Sample(z_n, v, t_n, sample_idx);
    h = 1e-6;
    J0 = md_cont_update.Discrepancy_Evaluation_Sample(z_n, t_n, sample_idx);
    J1 = md_cont_update.Discrepancy_Evaluation_Sample(z_n + h * v, t_n, sample_idx);
    Jv_fd = (J1 - J0) / (h);
    % fprintf('Rel. err in finite difference: %.3e\n', norm(Jv - Jv_fd)/norm(Jv))

    % Adjoint check for Apply_Discrepancy_z_Jacobian_transpose_Sample
    JT_w = md_cont_update.Apply_Discrepancy_z_Jacobian_transpose_Sample(z_n, w, t_n, sample_idx);
    lhs = Jv(:)' * w(:);
    rhs = v(:)' * JT_w(:);
    % fprintf('Rel. err in adjoint operator: %.3e\n', abs(lhs - rhs) / abs(lhs));

    % [u_k, z_k, beta_k] = md_cont_update.Posterior_Update_Sample(sample_idx);
end

% ---Checking if CG converges or not---
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
