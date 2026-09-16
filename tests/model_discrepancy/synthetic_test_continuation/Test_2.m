%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
clear;
close all;
rng(121234);

suppress_output = true;

m = 51;

data_interface = MD_Data_Interface_synthetic_test_continuation(m);
data_interface.Load_Data();

u_prior_interface = MD_u_Prior_Interface_synthetic_test_continuation(m);
z_prior_interface = MD_z_Prior_Interface_synthetic_test_continuation(m);

md_post_sampling = MD_Posterior_Sampling(data_interface, u_prior_interface, z_prior_interface);
alpha_d = 1.e-5;
num_post_samples = 10;
md_post_sampling.Compute_Posterior_Data(alpha_d, num_post_samples);

opt_prob_interface = MD_Opt_Prob_Interface_synthetic_test_continuation(m);
md_hessian_analysis = MD_Hessian_Analysis(opt_prob_interface, z_prior_interface);
num_evals = 30;
oversampling = 10;
md_hessian_analysis.Compute_Hessian_GEVP(data_interface.z_opt, num_evals, oversampling);

% Setup
sample_idx = 1;
r = length(data_interface.z_opt);
if ~isempty(md_hessian_analysis.evals)
    r = length(md_hessian_analysis.evals);
end

beta_n = randn(r, 1);
z_n = data_interface.z_opt + md_hessian_analysis.Apply_V(beta_n);
v_z = randn(size(data_interface.z_opt));
v_beta = randn(r, 1);
w_u = randn(size(data_interface.u_opt));
sen_op = MD_Continuation_Sensitivity_Operators(md_post_sampling, md_hessian_analysis);

error = [];
tol = 1.e-3;

% Check individual sample code
[delta_mean, ~] = md_post_sampling.Posterior_Discrepancy_Samples(z_n);
delta_mean = delta_mean{1};

delta_mean_1 = sen_op.Discrepancy_Evaluation_Mean(z_n);
delta_sample_1 = sen_op.Discrepancy_Evaluation_Sample_Beta(beta_n, sample_idx);

error(end + 1) = norm(delta_mean_1 - delta_mean) / max(1, norm(delta_mean));
if ~suppress_output
    fprintf('Rel. err in discrep @ mean:   %.3e\n', error(end));
end

sample_finite_err = double(~all(isfinite(delta_sample_1)));
error(end + 1) = sample_finite_err;
if ~suppress_output
    fprintf('Finite-output err for sample discrep: %.3e\n', error(end));
end

% Finite Difference check for Apply_Discrepancy_z_Jacobian_Mean
Jv = sen_op.Apply_Discrepancy_z_Jacobian_Mean(v_z);
h = 1e-6;
J0 = sen_op.Discrepancy_Evaluation_Mean(z_n);
J1 = sen_op.Discrepancy_Evaluation_Mean(z_n + h * v_z);
Jv_fd = (J1 - J0) / (h);
error(end + 1) = norm(Jv - Jv_fd) / max(1, norm(Jv));
if ~suppress_output
    fprintf('Rel. err in mean finite difference: %.3e\n', error(end));
end

% Adjoint check for Apply_Discrepancy_z_Jacobian_Transpose_Mean
JT_w = sen_op.Apply_Discrepancy_z_Jacobian_Transpose_Mean(w_u);
lhs = Jv(:)' * w_u(:);
rhs = v_z(:)' * JT_w(:);
error(end + 1) = abs(lhs - rhs) / max([1, abs(lhs), abs(rhs)]);
if ~suppress_output
    fprintf('Rel. err in mean adjoint operator: %.3e\n', error(end));
end

% Finite Difference check for Apply_Discrepancy_Beta_Jacobian_Sample
Jv = sen_op.Apply_Discrepancy_Beta_Jacobian_Sample(v_beta, sample_idx);
h = 1e-6;
J0 = sen_op.Discrepancy_Evaluation_Sample_Beta(beta_n, sample_idx);
J1 = sen_op.Discrepancy_Evaluation_Sample_Beta(beta_n + h * v_beta, sample_idx);
Jv_fd = (J1 - J0) / (h);
error(end + 1) = norm(Jv - Jv_fd) / max(1, norm(Jv));
if ~suppress_output
    fprintf('Rel. err in sample finite difference: %.3e\n', error(end));
end

% Adjoint check for Apply_Discrepancy_Beta_Jacobian_Transpose_Sample
JT_w = sen_op.Apply_Discrepancy_Beta_Jacobian_Transpose_Sample(w_u, sample_idx);
lhs = Jv(:)' * w_u(:);
rhs = v_beta(:)' * JT_w(:);
error(end + 1) = abs(lhs - rhs) / max([1, abs(lhs), abs(rhs)]);
if ~suppress_output
    fprintf('Rel. err in sample adjoint operator: %.3e\n', error(end));
end

if max(error) > tol
    fprintf(2, '\nmodel_discrepancy/synthetic_test_continuation Test_2 failed.\n');
else
    fprintf(1, '\nmodel_discrepancy/synthetic_test_continuation Test_2 passed.\n');
end
