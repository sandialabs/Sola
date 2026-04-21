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

% Setup
z_n = randn(size(data_interface.z_opt));
t_n = 1.0;
sample_idx = 1;
v = randn(size(data_interface.z_opt));
w = randn(size(data_interface.u_opt));
sen_op = MD_Continuation_Sensitivity_Operators(md_post_sampling, md_hessian_analysis);

% Check individual sample code
[delta_mean, delta_samples] = md_post_sampling.Posterior_Discrepancy_Samples(z_n);
delta_mean = delta_mean{1};
delta_samples = delta_samples{1};

delta_mean_1 = sen_op.Discrepancy_Evaluation_Mean(z_n);
delta_samples_1 = sen_op.Discrepancy_Evaluation_Sample(z_n, 1);
error = norm(delta_mean_1 - delta_mean) / norm(delta_mean);
if ~suppress_output
    fprintf('Rel. err in discrep @ mean:   %.3e\n', error(end));
end
error(end + 1) = norm(delta_samples_1 - delta_samples(:, 1)) / norm(delta_samples(:, 1));
if ~suppress_output
    fprintf('Rel. err in discrep @ sample:   %.3e\n', error(end));
end

% Finite Difference check for Apply_Discrepancy_z_Jacobian_Sample
Jv = sen_op.Apply_Discrepancy_z_Jacobian_Mean(v);
h = 1e-6;
J0 = sen_op.Discrepancy_Evaluation_Mean(z_n);
J1 = sen_op.Discrepancy_Evaluation_Mean(z_n + h * v);
Jv_fd = (J1 - J0) / (h);
error(end + 1) = norm(Jv - Jv_fd) / norm(Jv);
if ~suppress_output
    fprintf('Rel. err in finite difference: %.3e\n', error(end));
end

% Adjoint check for Apply_Discrepancy_z_Jacobian_transpose_Sample
JT_w = sen_op.Apply_Discrepancy_z_Jacobian_Transpose_Mean(w);
lhs = Jv(:)' * w(:);
rhs = v(:)' * JT_w(:);
error(end + 1) = abs(lhs - rhs) / abs(lhs);
if ~suppress_output
    fprintf('Rel. err in adjoint operator: %.3e\n', error(end));
end

% Finite Difference check for Apply_Discrepancy_z_Jacobian_Sample
Jv = sen_op.Apply_Discrepancy_z_Jacobian_Sample(v, z_n, sample_idx);
h = 1e-6;
J0 = sen_op.Discrepancy_Evaluation_Sample(z_n, sample_idx);
J1 = sen_op.Discrepancy_Evaluation_Sample(z_n + h * v, sample_idx);
Jv_fd = (J1 - J0) / (h);
error(end + 1) = norm(Jv - Jv_fd) / norm(Jv);
if ~suppress_output
    fprintf('Rel. err in finite difference: %.3e\n', error(end));
end

% Adjoint check for Apply_Discrepancy_z_Jacobian_transpose_Sample
JT_w = sen_op.Apply_Discrepancy_z_Jacobian_Transpose_Sample(w, z_n, sample_idx);
lhs = Jv(:)' * w(:);
rhs = v(:)' * JT_w(:);
error(end + 1) = abs(lhs - rhs) / abs(lhs);
if ~suppress_output
    fprintf('Rel. err in adjoint operator: %.3e\n', error(end));
end

if error > 1.e-3
    fprintf(2, '\nmodel_discrepancy/synthetic_test_continuation Test_2 failed.\n');
else
    fprintf(1, '\nmodel_discrepancy/synthetic_test_continuation Test_2 passed.\n');
end
