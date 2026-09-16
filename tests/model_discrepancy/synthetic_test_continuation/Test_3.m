clear;
close all;
rng(121235);

fprintf('\n============================================================\n');
fprintf('Running model_discrepancy/synthetic_test_continuation Test 3\n');
fprintf('============================================================\n\n');

%% Problem setup

m = 51;

data_interface = MD_Data_Interface_synthetic_test_continuation(m);
data_interface.Load_Data();

u_prior_interface = MD_u_Prior_Interface_synthetic_test_continuation(m);
z_prior_interface = MD_z_Prior_Interface_synthetic_test_continuation(m);

try
    wu_test_in = randn(m, 2);
    wu_test_out = u_prior_interface.Apply_W_u(wu_test_in);

    if isempty(wu_test_out) || any(size(wu_test_out) ~= size(wu_test_in))
        error('Apply_W_u returned an empty or incorrectly sized output.');
    end
catch ME
    fprintf(2, '\nTest setup failed: u_prior_interface.Apply_W_u is required.\n');
    fprintf(2, 'Add the following method to MD_u_Prior_Interface_synthetic_test_continuation:\n\n');
    fprintf(2, '    function [u_out] = Apply_W_u(this, u_in)\n');
    fprintf(2, '        u_out = this.W_u * u_in;\n');
    fprintf(2, '    end\n\n');
    rethrow(ME);
end

md_post_sampling = MD_Posterior_Sampling(data_interface, u_prior_interface, z_prior_interface);

alpha_d = 1.e-5;
num_post_samples = 5;
md_post_sampling.Compute_Posterior_Data(alpha_d, num_post_samples, false);

opt_prob_interface = MD_Opt_Prob_Interface_synthetic_test_continuation(m);

md_hessian_analysis = MD_Hessian_Analysis(opt_prob_interface, z_prior_interface);

num_evals = 10;
oversampling = 10;

md_hessian_analysis.Compute_Hessian_GEVP(data_interface.z_opt, num_evals, oversampling);

if isempty(md_hessian_analysis.evals)
    r = length(data_interface.z_opt);
else
    r = length(md_hessian_analysis.evals);
end

fprintf('Problem dimensions:\n');
fprintf('  m                 = %d\n', m);
fprintf('  reduced dim r     = %d\n', r);
fprintf('  posterior samples = %d\n\n', num_post_samples);

%% Construct continuation sensitivity operator

sen_op = MD_Continuation_Sensitivity_Operators(md_post_sampling, md_hessian_analysis);

sample_idx = 1;
t0 = 0.50;
time_index = 1;

theta_traj = make_theta_traj(t0, sample_idx);

beta = 5.e-2 * randn(r, 1);

v = randn(r, 1);
v = v / norm(v);

w = randn(r, 1);
w = w / norm(w);

u_test = randn(m, 1);

%% Test bookkeeping

test_names = {};
test_errs = [];
test_tols = [];

%% Lazy matrix-normal sampler validation

fprintf('Running lazy matrix-normal sampler validation tests...\n\n');

sampler_fwd_persist = fresh_breve_sampler(md_hessian_analysis, z_prior_interface, md_post_sampling.post_data,u_prior_interface, m);

v_persist = randn(r, 1);
u_enrich = randn(m, 1);

y_before = sampler_fwd_persist.Eval(v_persist);

sampler_fwd_persist.Apply_Jacobian_Transpose(u_enrich);

y_after = sampler_fwd_persist.Eval(v_persist);

forward_persistence_err = norm(y_after - y_before) / max(1, norm(y_before));

test_names{end+1} = 'Lazy sampler forward persistence after adjoint enrichment';
test_errs(end+1) = forward_persistence_err;
test_tols(end+1) = 1.e-11;

sampler_adj_persist = fresh_breve_sampler(md_hessian_analysis, z_prior_interface, md_post_sampling.post_data,u_prior_interface, m);

u_persist = randn(m, 1);
v_enrich = randn(r, 1);

z_before = sampler_adj_persist.Apply_Jacobian_Transpose(u_persist);

sampler_adj_persist.Eval(v_enrich);

z_after = sampler_adj_persist.Apply_Jacobian_Transpose(u_persist);

adjoint_persistence_err = norm(z_after - z_before) / max(1, norm(z_before));

test_names{end+1} = 'Lazy sampler adjoint persistence after forward enrichment';
test_errs(end+1) = adjoint_persistence_err;
test_tols(end+1) = 1.e-11;

sampler_interleave = fresh_breve_sampler(md_hessian_analysis, z_prior_interface, md_post_sampling.post_data,u_prior_interface, m);

v0_interleave = randn(r, 1);
u0_interleave = randn(m, 1);

y0_interleave = sampler_interleave.Eval(v0_interleave);
z0_interleave = sampler_interleave.Apply_Jacobian_Transpose(u0_interleave);

num_interleave_queries = 10;
for jj = 1:num_interleave_queries
    sampler_interleave.Eval(randn(r, 1));
    sampler_interleave.Apply_Jacobian_Transpose(randn(m, 1));
end

y1_interleave = sampler_interleave.Eval(v0_interleave);
z1_interleave = sampler_interleave.Apply_Jacobian_Transpose(u0_interleave);

interleave_forward_err = norm(y1_interleave - y0_interleave) / max(1, norm(y0_interleave));

interleave_adjoint_err = norm(z1_interleave - z0_interleave) / max(1, norm(z0_interleave));

test_names{end+1} = 'Lazy sampler interleaved forward repeatability';
test_errs(end+1) = interleave_forward_err;
test_tols(end+1) = 1.e-11;

test_names{end+1} = 'Lazy sampler interleaved adjoint repeatability';
test_errs(end+1) = interleave_adjoint_err;
test_tols(end+1) = 1.e-11;

sampler_explicit = fresh_breve_sampler(md_hessian_analysis, z_prior_interface, md_post_sampling.post_data,u_prior_interface, m);

fully_explore_breve_sampler(sampler_explicit);

X_cache = explicit_matrix_from_lazy_sampler(sampler_explicit);

v_explicit = randn(r, 1);
u_explicit = randn(m, 1);

lazy_forward = sampler_explicit.Eval(v_explicit);
explicit_forward = X_cache * v_explicit;

lazy_adjoint = sampler_explicit.Apply_Jacobian_Transpose(u_explicit);
explicit_adjoint = X_cache' * u_explicit;

explicit_forward_err = norm(lazy_forward - explicit_forward) / max(1, norm(explicit_forward));

explicit_adjoint_err = norm(lazy_adjoint - explicit_adjoint) / max(1, norm(explicit_adjoint));

test_names{end+1} = 'Lazy sampler reconstructed-matrix forward consistency';
test_errs(end+1) = explicit_forward_err;
test_tols(end+1) = 1.e-09;

test_names{end+1} = 'Lazy sampler reconstructed-matrix adjoint consistency';
test_errs(end+1) = explicit_adjoint_err;
test_tols(end+1) = 1.e-09;

num_moment_samples = 400;
num_functionals = 4;

A_moment = randn(m, num_functionals);
B_moment = randn(r, num_functionals);

phi_samples = zeros(num_moment_samples, num_functionals);

for ss = 1:num_moment_samples

    sampler_moment = fresh_breve_sampler(md_hessian_analysis, z_prior_interface, md_post_sampling.post_data, u_prior_interface, m);

    for jj = 1:num_functionals
        Xbj = sampler_moment.Eval(B_moment(:, jj));
        phi_samples(ss, jj) = A_moment(:, jj)' * Xbj;

        sampler_moment.Apply_Jacobian_Transpose(randn(m, 1));
    end

end

empirical_mean = mean(phi_samples, 1)';
empirical_cov = cov(phi_samples, 1);

predicted_cov = zeros(num_functionals, num_functionals);

sigma_ref_sampler = fresh_breve_sampler(md_hessian_analysis, z_prior_interface, md_post_sampling.post_data, u_prior_interface, m);

for jj = 1:num_functionals
    Wuinv_aj = u_prior_interface.Apply_W_u_Inverse(A_moment(:, jj));
    Sig_bj = sigma_ref_sampler.Apply_Sigma_Beta(B_moment(:, jj));

    for kk = 1:num_functionals
        left_cov = A_moment(:, kk)' * Wuinv_aj;
        right_cov = B_moment(:, kk)' * Sig_bj;
        predicted_cov(kk, jj) = left_cov * right_cov;
    end
end

mean_scale = sqrt(max(1, trace(predicted_cov)));
moment_mean_err = norm(empirical_mean) / mean_scale;

moment_cov_err = norm(empirical_cov - predicted_cov, 'fro') / max(1, norm(predicted_cov, 'fro'));

test_names{end+1} = 'Lazy sampler empirical scalar-functional mean';
test_errs(end+1) = moment_mean_err;
test_tols(end+1) = 2.5e-1;

test_names{end+1} = 'Lazy sampler empirical scalar-functional covariance';
test_errs(end+1) = moment_cov_err;
test_tols(end+1) = 3.5e-1;

%% Fully explore the breve sampler before derivative checks

sampler = sen_op.Get_Breve_Sampler(sample_idx);

fully_explore_breve_sampler(sampler);

[kl, kr] = sampler.lazy_X.Basis_Dimensions();

fprintf('Fully explored sample %d lazy breve sampler:\n', sample_idx);
fprintf('  left basis dim  = %d\n', kl);
fprintf('  right basis dim = %d\n\n', kr);

%% Breve sampler orthonormality and adjoint consistency

lazy_X = sampler.lazy_X;

if isempty(lazy_X.Q_l)
    err_left_orth = 0.0;
else
    WQl = u_prior_interface.Apply_W_u(lazy_X.Q_l);
    err_left_orth = norm(lazy_X.Q_l' * WQl - eye(size(lazy_X.Q_l, 2)), 'fro');
end

if isempty(lazy_X.Q_r)
    err_right_orth = 0.0;
else
    err_right_orth = norm(lazy_X.Q_r' * lazy_X.Sigma_Q_r - eye(size(lazy_X.Q_r, 2)), 'fro');
end

beta_adj = randn(r, 1);
u_adj = randn(m, 1);

Xbeta = sampler.Eval(beta_adj);
XTu = sampler.Apply_Jacobian_Transpose(u_adj);

breve_left = Xbeta' * u_adj;
breve_right = beta_adj' * XTu;

breve_adj_err = abs(breve_left - breve_right) / max([1, abs(breve_left), abs(breve_right)]);

test_names{end+1} = 'Lazy breve sampler left weighted orthonormality';
test_errs(end+1) = err_left_orth;
test_tols(end+1) = 1.e-10;

test_names{end+1} = 'Lazy breve sampler right Sigma_beta orthonormality';
test_errs(end+1) = err_right_orth;
test_tols(end+1) = 1.e-09;

test_names{end+1} = 'Lazy breve sampler adjoint relative error';
test_errs(end+1) = breve_adj_err;
test_tols(end+1) = 1.e-10;

if isempty(lazy_X.Q_l) || isempty(lazy_X.Q_r)
    cache_consistency_err = 0.0;
else
    left_cache = lazy_X.Q_l' * u_prior_interface.Apply_W_u(lazy_X.Y);
    right_cache = lazy_X.T' * lazy_X.Q_r;

    cache_consistency_err = norm(left_cache - right_cache, 'fro') / max(1, norm(right_cache, 'fro'));
end

test_names{end+1} = 'Lazy sampler revealed-action cache consistency';
test_errs(end+1) = cache_consistency_err;
test_tols(end+1) = 1.e-09;

%% Full beta-space sample discrepancy Jacobian finite difference

eps_fd_disc = 1.e-5;

D0 = sen_op.Discrepancy_Evaluation_Sample_Beta(beta, sample_idx);

D1 = sen_op.Discrepancy_Evaluation_Sample_Beta(beta + eps_fd_disc * v, sample_idx);

fd_D = (D1 - D0) / eps_fd_disc;

jac_D = sen_op.Apply_Discrepancy_Beta_Jacobian_Sample(v, sample_idx);

disc_jac_fd_err = norm(fd_D - jac_D) / max(1, norm(jac_D));

test_names{end+1} = 'Full beta sample discrepancy Jacobian FD relative error';
test_errs(end+1) = disc_jac_fd_err;
test_tols(end+1) = 1.e-5;

%% Full beta-space sample discrepancy adjoint consistency

Dv = sen_op.Apply_Discrepancy_Beta_Jacobian_Sample(v, sample_idx);
DTu = sen_op.Apply_Discrepancy_Beta_Jacobian_Transpose_Sample(u_test, sample_idx);

disc_left = Dv' * u_test;
disc_right = v' * DTu;

disc_adj_err = abs(disc_left - disc_right) / max([1, abs(disc_left), abs(disc_right)]);

test_names{end+1} = 'Full beta sample discrepancy adjoint relative error';
test_errs(end+1) = disc_adj_err;
test_tols(end+1) = 1.e-10;

%% Sample continuation gradient finite difference

eps_fd_grad = 1.e-6;

[g0, val0] = sen_op.Gradient(beta, theta_traj, time_index);

[~, val1] = sen_op.Gradient(beta + eps_fd_grad * v, theta_traj, time_index);

fd_grad = (val1 - val0) / eps_fd_grad;
dir_grad = g0' * v;

grad_fd_err = abs(fd_grad - dir_grad) / max([1, abs(dir_grad), abs(fd_grad)]);

test_names{end+1} = 'Sample continuation gradient FD relative error';
test_errs(end+1) = grad_fd_err;
test_tols(end+1) = 5.e-4;

%% Sample continuation Hessian finite difference

eps_fd_hess = 1.e-6;

[g0, ~] = sen_op.Gradient(beta, theta_traj, time_index);

[g1, ~] = sen_op.Gradient(beta + eps_fd_hess * v, theta_traj, time_index);

fd_Hv = (g1 - g0) / eps_fd_hess;

Hv = sen_op.Apply_Hessian(v, beta, theta_traj, time_index);

hess_fd_err = norm(fd_Hv - Hv) / max(1, norm(Hv));

test_names{end+1} = 'Sample continuation Hessian FD relative error';
test_errs(end+1) = hess_fd_err;
test_tols(end+1) = 5.e-3;

%% Sample continuation Hessian symmetry

Hv = sen_op.Apply_Hessian(v, beta, theta_traj, time_index);
Hw = sen_op.Apply_Hessian(w, beta, theta_traj, time_index);

hess_left = v' * Hw;
hess_right = w' * Hv;

hess_sym_err = abs(hess_left - hess_right) / max([1, abs(hess_left), abs(hess_right)]);

test_names{end+1} = 'Sample continuation Hessian symmetry relative error';
test_errs(end+1) = hess_sym_err;
test_tols(end+1) = 1.e-7;

%% Sample continuation mixed derivative Apply_B finite difference

eps_fd_B = 1.e-6;

theta_traj_t0 = make_theta_traj(t0, sample_idx);
theta_traj_t1 = make_theta_traj(t0 + eps_fd_B, sample_idx);

[g_t0, ~] = sen_op.Gradient(beta, theta_traj_t0, time_index);
[g_t1, ~] = sen_op.Gradient(beta, theta_traj_t1, time_index);

fd_B = (g_t1 - g_t0) / eps_fd_B;

B = sen_op.Apply_B(beta, theta_traj_t0, time_index);

B_fd_err = norm(fd_B - B) / max(1, norm(B));

test_names{end+1} = 'Sample continuation Apply_B FD relative error';
test_errs(end+1) = B_fd_err;
test_tols(end+1) = 5.e-4;

%% End-to-end posterior sample continuation finite-output check

num_continuation_steps = 2;

md_cont_update = MD_Continuation_Update(md_post_sampling, md_hessian_analysis, num_continuation_steps);

[u_ks, z_ks, beta_ks] = md_cont_update.Posterior_Update_Samples();

finite_ok = all(isfinite(u_ks(:))) && all(isfinite(z_ks(:))) && all(isfinite(beta_ks(:)));

if finite_ok
    finite_err = 0.0;
else
    finite_err = Inf;
end

test_names{end+1} = 'End-to-end posterior sample continuation finite outputs';
test_errs(end+1) = finite_err;
test_tols(end+1) = 0.0;

%% Report results

fprintf('\nTest results:\n');
fprintf('------------------------------------------------------------\n');

all_passed = true;

for j = 1:length(test_names)
    passed = test_errs(j) <= test_tols(j);

    if passed
        status = 'PASS';
    else
        status = 'FAIL';
        all_passed = false;
    end

    fprintf('%-65s  err = %.4e   tol = %.4e   %s\n',  test_names{j}, test_errs(j), test_tols(j), status);
end

fprintf('------------------------------------------------------------\n');

if all_passed
    fprintf(1, '\nmodel_discrepancy/synthetic_test_continuation 3 passed.\n\n');
else
    fprintf(2, '\nmodel_discrepancy/synthetic_test_continuation Test 3 failed.\n\n');
end

%% Local helpers

function theta_traj = make_theta_traj(t, sample_idx)

    theta_traj.Get_Time = @(time_index) t;
    theta_traj.Get_Sample_Index = @() sample_idx;

end

function fully_explore_breve_sampler(sampler)

    lazy_X = sampler.lazy_X;

    input_dim = lazy_X.input_dim;
    output_dim = lazy_X.output_dim;

    for j = 1:input_dim
        e = zeros(input_dim, 1);
        e(j) = 1.0;
        lazy_X.Forward_Apply(e);
    end

    for i = 1:output_dim
        e = zeros(output_dim, 1);
        e(i) = 1.0;
        lazy_X.Adjoint_Apply(e);
    end

end

function sampler = fresh_breve_sampler(md_hessian_analysis, z_prior_interface, post_data, u_prior_interface, output_dim)

    lazy_tol = 1.e-12;

    sampler = MD_Breve_Beta_Sampler(md_hessian_analysis, z_prior_interface, post_data, u_prior_interface, output_dim, lazy_tol);

end

function X_cache = explicit_matrix_from_lazy_sampler(sampler)

    lazy_X = sampler.lazy_X;

    input_dim = lazy_X.input_dim;

    Sigma_beta = zeros(input_dim, input_dim);

    for j = 1:input_dim
        e = zeros(input_dim, 1);
        e(j) = 1.0;
        Sigma_beta(:, j) = sampler.Apply_Sigma_Beta(e);
    end

    if isempty(lazy_X.Q_r)
        X_cache = zeros(lazy_X.output_dim, lazy_X.input_dim);
        return;
    end

    X_cache = lazy_X.Y * (lazy_X.Q_r' * Sigma_beta);

end