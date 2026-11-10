%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Test of posterior sample completion (paper: "Posterior sample completion")
% for the lazy breve sampler X ~ MN(0, W_u^{-1}, Sigma_beta).
%
%   Test 1: completion reproduces the revealed data and keeps forward and
%           adjoint applies mutually consistent, before and after completion.
%   Test 2: completed samples have the same law as lazy samples (Monte Carlo
%           comparison of E||X b||^2 for a fresh direction b, and of the
%           conditional mean of X b given the revealed data).
%   Test 3: automatic switch to the completed sample at a cache-size threshold.
%   Test 4: end-to-end continuation sampling with and without completion.
%   Test 5: measured storage and runtime of lazy, completed, and explicit sampling.
%   Test 6: same comparison inside actual mean-to-sample continuation.

clearvars;
addpath(genpath('..'));
OED_Setup;
oed_results = load('oed-results.mat');
num_data_points = 3;
data_interface.Set_Z_and_D(oed_results.Z_oed(:, 1:num_data_points), oed_results.D_oed(:, 1:num_data_points));

% ------------------------
% User-adjustable settings
% ------------------------
num_samples = 60;       % posterior samples (also used in Test 4)
num_prior_queries = 5;      % distinct queries revealed before completion
num_mc = 400;               % Monte Carlo repetitions for Test 2
run_end_to_end = false;      % Test 4 (runs continuation)
pass_tol = 1e-6;
rng(2468);

md_post_sampling = MD_Posterior_Sampling(data_interface, u_prior_interface, z_prior_interface);
md_post_sampling.Compute_Posterior_Data(alpha_d, num_samples, true);
sen_op = MD_Continuation_Sensitivity_Operators(md_post_sampling, md_hessian_analysis, false, false, 0);
sampler = sen_op.Get_Breve_Sampler(1);
n_u = sampler.output_dim;
r = sampler.input_dim;
fprintf('\nn_u = %d, r = %d\n', n_u, r);

results = struct();

%% Test 1: reproduction of revealed data and consistency
fprintf('\n--- Test 1: data reproduction and consistency ---\n');
lazy = sampler.lazy_X;
B = randn(r, num_prior_queries);
U = u_prior_interface.Sample_with_Covariance_W_u_Inverse(num_prior_queries);
for j = 1:num_prior_queries
    lazy.Forward_Apply(B(:, j));
    lazy.Adjoint_Apply(U(:, j));
end
[kl, kr] = lazy.Basis_Dimensions();
fprintf('Cache sizes before completion: k_f = %d, k_t = %d\n', kr, kl);

Y_before = zeros(n_u, num_prior_queries);
X_T_before = zeros(r, num_prior_queries);
for j = 1:num_prior_queries
    Y_before(:, j) = lazy.Forward_Apply(B(:, j));
    X_T_before(:, j) = lazy.Adjoint_Apply(U(:, j));
end

lazy.Complete_Sample();
assert(lazy.is_completed);
fprintf('S has %d columns (rank of residual right covariance)\n', size(lazy.S_right, 2));

Y_after = zeros(n_u, num_prior_queries);
X_T_after = zeros(r, num_prior_queries);
for j = 1:num_prior_queries
    Y_after(:, j) = lazy.Forward_Apply(B(:, j));
    X_T_after(:, j) = lazy.Adjoint_Apply(U(:, j));
end
err_fwd = norm(Y_after - Y_before, 'fro') / max(norm(Y_before, 'fro'), eps);
err_adj = norm(X_T_after - X_T_before, 'fro') / max(norm(X_T_before, 'fro'), eps);

% Cache identities X Q_r = Y and Q_l' W_u X = T'
XQr = zeros(n_u, kr);
for j = 1:kr
    XQr(:, j) = lazy.Forward_Apply(lazy.Q_r(:, j));
end
err_Y = norm(XQr - lazy.Y, 'fro') / max(norm(lazy.Y, 'fro'), eps);
err_T = 0;
for j = 1:kl
    t_j = lazy.Adjoint_Apply(lazy.W_Q_l(:, j));
    err_T = max(err_T, norm(t_j - lazy.T(:, j)) / max(norm(lazy.T(:, j)), eps));
end

% Mutual consistency on fresh queries, and against a materialized matrix
b_new = randn(r, 1);
u_new = u_prior_interface.Sample_with_Covariance_W_u_Inverse(1);
Xb = lazy.Forward_Apply(b_new);
Xtu = lazy.Adjoint_Apply(u_new);
err_cons = abs(u_new' * Xb - Xtu' * b_new) / max(abs(u_new' * Xb), eps);
X_full = zeros(n_u, r);
for j = 1:r
    e_j = zeros(r, 1); e_j(j) = 1;
    X_full(:, j) = lazy.Forward_Apply(e_j);
end
err_mat = max(norm(X_full * b_new - Xb) / norm(Xb), norm(X_full' * u_new - Xtu) / norm(Xtu));
fprintf('Rel. change in earlier forward queries : %.2e\n', err_fwd);
fprintf('Rel. change in earlier adjoint queries : %.2e\n', err_adj);
fprintf('Cache identity X Q_r = Y               : %.2e\n', err_Y);
fprintf('Cache identity Q_l'' W_u X = T''         : %.2e\n', err_T);
fprintf('Fresh forward/adjoint consistency      : %.2e\n', err_cons);
fprintf('Matrix-free vs materialized X          : %.2e\n', err_mat);
results.test1_pass = max([err_fwd, err_adj, err_Y, err_T, err_cons, err_mat]) < pass_tol;
fprintf('Test 1: %s\n', pass_str(results.test1_pass));

%% Test 2: statistical equivalence of lazy and completed samples
fprintf('\n--- Test 2: law of completed vs lazy samples (%d trials) ---\n', num_mc);
b_probe = randn(r, 1);
B_rev = randn(r, num_prior_queries);
val_lazy = zeros(num_mc, 1);
val_comp = zeros(num_mc, 1);
rev_comp = zeros(num_mc, 1);
for k = 1:num_mc
    s_lazy = make_fresh_lazy(sen_op, sampler);
    s_comp = make_fresh_lazy(sen_op, sampler);
    for j = 1:num_prior_queries
        s_lazy.Forward_Apply(B_rev(:, j));
        s_comp.Forward_Apply(B_rev(:, j));
    end
    rev_before = s_comp.Forward_Apply(B_rev(:, 1));
    s_comp.Complete_Sample();
    rev_comp(k) = norm(s_comp.Forward_Apply(B_rev(:, 1)) - rev_before) / norm(rev_before);
    val_lazy(k) = norm(s_lazy.Forward_Apply(b_probe))^2;
    val_comp(k) = norm(s_comp.Forward_Apply(b_probe))^2;
end
m_l = mean(val_lazy); m_c = mean(val_comp);
se = sqrt(var(val_lazy) / num_mc + var(val_comp) / num_mc);
z_stat = (m_c - m_l) / se;
fprintf('E||Xb||^2 lazy = %.4e, completed = %.4e (z = %.2f)\n', m_l, m_c, z_stat);
fprintf('Max change in revealed query on completion: %.2e\n', max(rev_comp));
results.test2_pass = abs(z_stat) < 4 && max(rev_comp) < pass_tol;
fprintf('Test 2: %s\n', pass_str(results.test2_pass));

%% Test 3: automatic switching at a cache-size threshold
fprintf('\n--- Test 3: automatic completion ---\n');
s_auto = make_fresh_lazy(sen_op, sampler);
s_auto.auto_complete = 0.05;   % complete at 5%% of break-even (forces a switch)
nq = 0;
while ~s_auto.is_completed && nq < 4 * r
    s_auto.Forward_Apply(randn(r, 1));
    s_auto.Adjoint_Apply(u_prior_interface.Sample_with_Covariance_W_u_Inverse(1));
    nq = nq + 1;
end
[kl, kr] = s_auto.Basis_Dimensions();
fprintf('Switched after %d query pairs with k_f = %d, k_t = %d\n', nq, kr, kl);
sizes_at_switch = [kl, kr];
for j = 1:10
    s_auto.Forward_Apply(randn(r, 1));
    s_auto.Adjoint_Apply(u_prior_interface.Sample_with_Covariance_W_u_Inverse(1));
end
[kl, kr] = s_auto.Basis_Dimensions();
results.test3_pass = s_auto.is_completed && isequal([kl, kr], sizes_at_switch);
fprintf('Caches frozen after switch: %d\n', isequal([kl, kr], sizes_at_switch));
fprintf('Test 3: %s\n', pass_str(results.test3_pass));

%% Test 4: end-to-end continuation with and without completion
if run_end_to_end
    fprintf('\n--- Test 4: continuation sampling with/without completion ---\n');
    md_cont_update = MD_Continuation_Update(md_post_sampling, md_hessian_analysis);
    mean_options = struct();
    mean_options.num_steps = 7;
    mean_options.num_max_newton_steps = 5;
    mean_options.newton_gradient_tol = 1e-8;
    [u_mean, z_mean, beta_mean] = md_cont_update.Posterior_Update_Mean(mean_options);

    sample_options = struct();
    sample_options.num_mean_to_sample_steps = 3;
    sample_options.num_max_newton_steps = 5;
    sample_options.newton_gradient_tol = 1e-8;
    sample_options.use_parallel_sampling = false;
    sample_options.materialize_breve_samples = false;

    sample_options.complete_breve_samples = 0;
    [~, z_lazy, ~] = md_cont_update.Posterior_Update_Samples_From_Mean(u_mean, z_mean, beta_mean, sample_options);
    sample_options.complete_breve_samples = 0.05;
    [~, z_comp, ~] = md_cont_update.Posterior_Update_Samples_From_Mean(u_mean, z_mean, beta_mean, sample_options);

    sd_lazy = median(std(z_lazy, 0, 2));
    sd_comp = median(std(z_comp, 0, 2));
    fprintf('Median coordinate std  lazy = %.4e, completed = %.4e (ratio %.2f)\n', sd_lazy, sd_comp, sd_comp / sd_lazy);
    % Both runs share the same posterior data, so samples are comparable pairwise.
    pair_diff = norm(z_lazy - z_comp, 'fro') / norm(z_lazy - z_mean, 'fro');
    fprintf('Pairwise ||z_lazy - z_comp|| / ||z_lazy - z_mean|| = %.3f\n', pair_diff);
    fprintf('Mean-shift (rel) lazy = %.3e, completed = %.3e\n', ...
        norm(mean(z_lazy, 2) - z_mean) / norm(z_mean), norm(mean(z_comp, 2) - z_mean) / norm(z_mean));
    results.test4_pass = all(isfinite(z_comp(:))) && sd_comp / sd_lazy > 0.5 && sd_comp / sd_lazy < 2;
    fprintf('Test 4: %s (statistical sanity check only; samples differ in randomness)\n', pass_str(results.test4_pass));
end

%% Test 5: storage and runtime of lazy vs completed vs explicit sampling
fprintf('\n--- Test 5: storage / runtime benchmark ---\n');
num_query_pairs = 100;     % forward+adjoint pairs per sampler (random, i.e. independent directions)
num_reps = 20;
Bq = randn(r, num_query_pairs);
Uq = u_prior_interface.Sample_with_Covariance_W_u_Inverse(num_query_pairs);
modes = {'lazy', 'auto_complete (f=1, break-even)', 'explicit (materialized)'};
bench = struct('mode', modes, 'storage', 0, 'time', 0, 'switch_pairs', NaN);
for m = 1:numel(modes)
    t_acc = zeros(num_reps, 1);
    for rep = 1:num_reps
        fresh = MD_Breve_Beta_Sampler(sampler.hessian_analysis, sampler.z_prior_interface, ...
            sampler.post_data, sen_op.post_sampling.u_prior_interface, n_u, ...
            sen_op.lazy_sampling_tol, m == 3, double(m == 2));
        tic;
        for j = 1:num_query_pairs
            fresh.Apply_Jacobian(Bq(:, j));
            fresh.Apply_Jacobian_Transpose(Uq(:, j));
            if m == 2 && isnan(bench(m).switch_pairs) && fresh.lazy_X.is_completed
                bench(m).switch_pairs = j;
            end
        end
        t_acc(rep) = toc;
    end
    bench(m).time = median(t_acc);
    if m == 3
        bench(m).storage = numel(fresh.materialized_X);
    else
        lx = fresh.lazy_X;
        bench(m).storage = numel(lx.Q_l) + numel(lx.W_Q_l) + numel(lx.Q_r) + numel(lx.Sigma_Q_r) + ...
            numel(lx.Y) + numel(lx.T) + numel(lx.S_left) + numel(lx.S_right);
    end
    fprintf('%-34s storage = %6d doubles (%.2fx explicit), time = %.4f s per %d pairs\n', ...
        bench(m).mode, bench(m).storage, bench(m).storage / (n_u * r), bench(m).time, num_query_pairs);
end
fprintf('Auto-complete switched after pair %d\n', bench(2).switch_pairs);
results.test5_info = bench;

%% Test 6: storage and runtime on an actual continuation problem
fprintf('\n--- Test 6: continuation benchmark (lazy vs completed vs explicit) ---\n');
num_bench_samples = min(num_samples, 20);
num_mean_to_sample_steps = 3;
if ~exist('beta_mean', 'var')
    md_cont_update = MD_Continuation_Update(md_post_sampling, md_hessian_analysis);
    mean_options = struct();
    mean_options.num_steps = 7;
    mean_options.num_max_newton_steps = 5;
    mean_options.newton_gradient_tol = 1e-8;
    [~, ~, beta_mean] = md_cont_update.Posterior_Update_Mean(mean_options);
end
md_cont_update.options.samples.num_mean_to_sample_steps = num_mean_to_sample_steps;
md_cont_update.options.samples.num_max_newton_steps = 5;
md_cont_update.options.samples.newton_gradient_tol = 1e-8;
md_cont_update.options.samples.time_step.use = false;

modes6 = {'lazy', 'lazy + completion (break-even)', 'explicit (materialized)'};
bench6 = struct('mode', modes6, 'total_storage', 0, 'max_storage', 0, 'time', 0, ...
    'median_grad_norm', 0, 'median_kf_kt', [0 0], 'num_completed', 0); % median_kf_kt = [k_t, k_f]
for m = 1:numel(modes6)
    % discard_cache = false so that all samplers stay alive and can be measured
    sen_b = MD_Continuation_Sensitivity_Operators(md_post_sampling, md_hessian_analysis, false, m == 3, double(m == 2));
    grad_norms6 = zeros(num_bench_samples, 1);
    store6 = zeros(num_bench_samples, 1);
    kk6 = zeros(num_bench_samples, 2);
    done6 = false(num_bench_samples, 1);
    tic;
    for i = 1:num_bench_samples
        [~, ~, ~, ~, ~, ~, grad_norms6(i), ~] = md_cont_update.Posterior_Update_Sample_From_Mean_Core(beta_mean, i, sen_b);
    end
    bench6(m).time = toc;
    for i = 1:num_bench_samples
        smp = sen_b.Get_Breve_Sampler(i);
        store6(i) = sampler_storage(smp);
        [kk6(i, 1), kk6(i, 2)] = smp.lazy_X.Basis_Dimensions();
        done6(i) = smp.lazy_X.is_completed;
    end
    bench6(m).total_storage = sum(store6);
    bench6(m).max_storage = max(store6);
    bench6(m).median_grad_norm = median(grad_norms6);
    bench6(m).median_kf_kt = median(kk6, 1);
    bench6(m).num_completed = sum(done6);
    fprintf('%-32s storage/sample = %6d doubles (%.2fx explicit), total time = %.2f s (%d samples), median ||grad|| = %.2e, median (k_t,k_f) = (%g,%g), completed: %d\n', ...
        bench6(m).mode, round(bench6(m).total_storage / num_bench_samples), ...
        bench6(m).total_storage / num_bench_samples / (n_u * r), bench6(m).time, num_bench_samples, ...
        bench6(m).median_grad_norm, bench6(m).median_kf_kt(1), bench6(m).median_kf_kt(2), bench6(m).num_completed);
end
results.test6_info = bench6;

fprintf('\n=== Summary ===\n');
disp(results);

%% Helpers
function n = sampler_storage(smp)
    % Number of doubles stored by a breve sampler (caches, completion factors, or explicit X)
    if smp.Is_Materialized()
        n = numel(smp.materialized_X);
        return;
    end
    lx = smp.lazy_X;
    n = numel(lx.Q_l) + numel(lx.W_Q_l) + numel(lx.Q_r) + numel(lx.Sigma_Q_r) + ...
        numel(lx.Y) + numel(lx.T) + numel(lx.S_left) + numel(lx.S_right);
end

function s = pass_str(tf)
    if tf, s = 'PASS'; else, s = 'FAIL'; end
end

function lazy = make_fresh_lazy(sen_op, ref_sampler)
    % Fresh, independent lazy sampler with the same covariances as ref_sampler.
    fresh = MD_Breve_Beta_Sampler(ref_sampler.hessian_analysis, ref_sampler.z_prior_interface, ...
        ref_sampler.post_data, sen_op.post_sampling.u_prior_interface, ref_sampler.output_dim, ...
        sen_op.lazy_sampling_tol, false, 0);
    lazy = fresh.lazy_X;
end
