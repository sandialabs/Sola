%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Experiment 2c: is the curvature limit (eps_curv) of the a priori time-step
% heuristic needed, and how should it be set?
%
% Part 1 (ablation).  Mean-to-sample continuation is run with
%   - 'backtrack only': no a priori step size; start with a full step and shrink only
%     when the backtracking safeguard (gradient reduction, negative curvature) rejects it,
%   - fixed N_c steps (the baseline),
%   - 'both'      : eps_nl = eps_nl_both,  eps_curv swept (Inf = nonlinearity limit only),
%   - 'curv only' : eps_nl = Inf,          eps_curv swept (curvature limit only).
% Every run ends with Newton polishing that stops, per sample, once
% ||g_beta|| <= newton_gradient_tol (at most num_max_newton_steps steps), so the
% endpoint gradient is small by construction.  What differs is
%   cost = continuation steps (+ backtracks) + newton_cost * Newton steps taken
% and whether the sample reaches the tolerance at all ("converged"; a sample
% whose continuation ends outside Newton's basin does not).  Adaptive settings are
% compared with the best uniform fixed N_c chosen in hindsight and with the
% per-sample oracle (cheapest fixed N_c for each sample).
%
% Part 2 (predictor-error test).  For each sample, one Euler predictor step
% from the mean start at dt = 0.5 and 0.25.  Its gradient defect
% ||g(pred) - g(start)|| should scale like dt^2; the coefficient c is compared
% (Spearman) with eta, eta * ||beta_dot||, A, A * ||beta_dot|| and ||beta_dot||,
% and with the difficulty of the sample (Newton steps needed after one
% continuation step).  If eta * ||beta_dot|| predicts c better than eta, the
% curvature limit lacks a path-length factor.

ci_num_data_points = 1;                % this experiment uses 1 data point (other drivers use 3)
CI_Setup;
CI_Sample_Reference;
clear ci_num_data_points;

fixed_sweep = [1, 2, 3, 4, 5];
num_max_newton_steps = 10;             % Newton polishing at t = 1: at most this many steps ...
newton_gradient_tol = lowfi_grad_norm; % ... stopping per sample once ||g_beta|| <= this
newton_cost = 0.5;                     % cost of a Newton step relative to a continuation step
eps_nl_both = 0.01;                    % eps_nl for the 'both' family
eps_curv_sweep = [1, 3, 10, 30, Inf];  % Inf = curvature limit off
pred_dts = [0.5, 0.25];

% Adaptive configurations.
configs = struct('family', {}, 'eps_nl', {}, 'eps_curv', {});
for ec = eps_curv_sweep
    configs(end + 1) = struct('family', 'both', 'eps_nl', eps_nl_both, 'eps_curv', ec); %#ok<SAGROW>
end
for ec = eps_curv_sweep(isfinite(eps_curv_sweep))
    configs(end + 1) = struct('family', 'curv only', 'eps_nl', Inf, 'eps_curv', ec); %#ok<SAGROW>
end
configs(end + 1) = struct('family', 'backtrack only', 'eps_nl', Inf, 'eps_curv', Inf);   % no a priori step: dt = 1, shrink on rejection
num_cfg = length(configs);
num_N = length(fixed_sweep);
S = num_samples;

fprintf('\n============================================================\n');
fprintf('Experiment 2c: curvature-limit sweep (%d data points, %d samples)\n', num_data_points, S);
fprintf('============================================================\n');
fprintf('Newton polishing: at most %d steps, stop at ||g_beta|| <= %.3e (low-fidelity optimizer level)\n', num_max_newton_steps, newton_gradient_tol);
fprintf('cost = continuation steps (+ backtracks) + %.2f * Newton steps; "converged" = reached the tolerance.\n', newton_cost);

%% Fixed-step baseline
fixed = cell(1, num_N);
fprintf('\nFixed steps (mean-to-sample):\n');
fprintf('%-28s | cont mean/max | Newton mean/max | cost mean / p90 / max | converged | s/sample\n', 'setting');
for jj = 1:num_N
    run_out = Run_Samples(sen_op, md_hessian_analysis, z_opt, r, 'fixed_mts', fixed_sweep(jj), beta_bar, [], [], sample_indices, num_max_newton_steps, newton_gradient_tol);
    run_out.summary = Summarize_CI_Cost(run_out, newton_cost, newton_gradient_tol, sprintf('fixed N_c = %d', fixed_sweep(jj)));
    fixed{jj} = run_out;
    Print_Summary(run_out.summary);
end
fixed_cost = cell2mat(cellfun(@(f) f.summary.cost, fixed, 'UniformOutput', false));   % S x num_N
fixed_converged = cellfun(@(f) f.summary.converged_frac, fixed);
fixed_mean_cost = cellfun(@(f) f.summary.mean_cost, fixed);
fixed_p90_cost = cellfun(@(f) f.summary.p90_cost, fixed);
ok_fixed = find(fixed_converged >= max(fixed_converged) - 1.e-12);
[best_fixed_cost, jb] = min(fixed_mean_cost(ok_fixed));
jb = ok_fixed(jb);
oracle_cost = mean(min(fixed_cost, [], 2), 'omitnan');
[~, jopt] = min(fixed_cost, [], 2);
N_opt = fixed_sweep(jopt)';            % cheapest fixed N_c for each sample (hindsight)
fprintf('Best uniform fixed N_c in hindsight (among those converging most often): N_c = %d, mean cost %.2f, p90 %.2f\n', ...
        fixed_sweep(jb), best_fixed_cost, fixed_p90_cost(jb));
fprintf('Per-sample oracle (cheapest fixed N_c per sample): mean cost %.2f\n', oracle_cost);

%% Adaptive settings
adapt = cell(1, num_cfg);
fprintf('\nAdaptive (Inf = limit off):\n');
fprintf('%-28s | cont mean/max | Newton mean/max | cost mean / p90 / max | converged | s/sample | curv binds | mean cost vs best fixed | p90 vs best fixed\n', 'setting');
for cc = 1:num_cfg
    cfg = configs(cc);
    bt_options = sample_bt_options;
    bt_options.epsilon_curv = cfg.eps_curv;
    if strcmp(cfg.family, 'backtrack only')
        run_out = Run_Samples(sen_op, md_hessian_analysis, z_opt, r, 'backtrack_mts', [], beta_bar, [], bt_options, sample_indices, num_max_newton_steps, newton_gradient_tol);
        label = 'backtrack only (dt = 1 start)';
    else
        run_out = Run_Samples(sen_op, md_hessian_analysis, z_opt, r, 'adaptive_mts', [], beta_bar, cfg.eps_nl, bt_options, sample_indices, num_max_newton_steps, newton_gradient_tol);
        label = sprintf('%s: nl=%g, curv=%g', cfg.family, cfg.eps_nl, cfg.eps_curv);
    end
    run_out.summary = Summarize_CI_Cost(run_out, newton_cost, newton_gradient_tol, label);
    run_out.family = cfg.family;  run_out.eps_nl = cfg.eps_nl;  run_out.eps_curv = cfg.eps_curv;
    run_out.curv_binds = mean(run_out.dt_curv < run_out.dt_nl, 'omitnan');   % fraction limited by curvature
    if strcmp(cfg.family, 'backtrack only')
        run_out.curv_binds = NaN;
    end
    adapt{cc} = run_out;
    sm = run_out.summary;
    fprintf('%-28s |  %4.2f / %2d   |   %4.2f / %2d     |   %5.2f / %5.2f / %5.2f |   %5.1f%%  | %.2f     |  %5.1f%%    |  %5.2fx                 |  %5.2fx\n', ...
            sm.label, sm.mean_cont, max(sm.cont_steps), sm.mean_newton, max(sm.newton_steps), ...
            sm.mean_cost, sm.p90_cost, sm.max_cost, 100 * sm.converged_frac, sm.time_per_sample, ...
            100 * run_out.curv_binds, sm.mean_cost / best_fixed_cost, sm.p90_cost / fixed_p90_cost(jb));
end
fprintf('(ratios < 1: cheaper than the best uniform fixed N_c chosen in hindsight; the oracle mean cost %.2f is the floor for any fixed-N_c rule.)\n', oracle_cost);

%% Part 2: predictor-error test
fprintf('\n------------------------------------------------------------\n');
fprintf('Predictor-error test: ||g(pred) - g(start)|| at dt = %s (expected ~ c dt^2)\n', mat2str(pred_dts));
fprintf('------------------------------------------------------------\n');
pred = struct('speed', nan(S, 1), 'c', nan(S, length(pred_dts)), 'slope', nan(S, 1));
for ii = 1:S
    idx = sample_indices(ii);
    try
        qn_prec = MD_Quasi_Newton_Preconditioner(md_hessian_analysis);
        qn_prec.Set_N(1);
        pt = Pseudo_Time_Continuation(beta_bar, sen_op, qn_prec);
        theta = MD_Discrepancy_Mean_To_Sample_Trajectory(1, idx);
        g0 = sen_op.Gradient(beta_bar, theta, 0);
        B = sen_op.Apply_B(beta_bar, theta, 0);
        f = pt.Apply_Inverse_Hessian(B, beta_bar, theta, 0);
        pred.speed(ii) = norm(f);
        defect = zeros(1, length(pred_dts));
        for kk = 1:length(pred_dts)
            g_pred = sen_op.Gradient(beta_bar - pred_dts(kk) * f, theta, pred_dts(kk));
            defect(kk) = norm(g_pred - g0);
            pred.c(ii, kk) = defect(kk) / pred_dts(kk)^2;
        end
        pred.slope(ii) = log(defect(1) / defect(end)) / log(pred_dts(1) / pred_dts(end));
    catch ME
        fprintf('  predictor test failed for sample %d: %s\n', idx, ME.message);
    end
end

j1 = find(fixed_sweep == 1, 1);
eta_s = adapt{1}.eta;                  % a priori quantities do not depend on eps
A_s = adapt{1}.A;
if isempty(j1)
    difficulty = nan(S, 1);
else
    difficulty = fixed{j1}.summary.newton_steps;   % Newton steps needed after one continuation step
end
c_pred = pred.c(:, end);

fprintf('\nPer-sample table:\n');
fprintf(' sample |    A       |   eta      | ||beta_dot|| | c (pred)   | slope | Newton after N_c=1 | cheapest N_c\n');
for ii = 1:S
    fprintf(' %5d  | %.3e | %.3e | %.3e   | %.3e | %5.2f |        %3d         |     %d\n', ii, A_s(ii), eta_s(ii), ...
            pred.speed(ii), c_pred(ii), pred.slope(ii), difficulty(ii), N_opt(ii));
end
fprintf('(slope should be ~2 for a clean quadratic predictor defect; cheapest N_c = argmin of the fixed-step cost per sample.)\n');

fprintf('\nSpearman correlation with the predictor defect coefficient c and with the sample difficulty:\n');
cands = {'eta', eta_s; 'eta*||beta_dot||', eta_s .* pred.speed; 'A', A_s; 'A*||beta_dot||', A_s .* pred.speed; '||beta_dot||', pred.speed};
pred_corr = zeros(size(cands, 1), 2);
for kk = 1:size(cands, 1)
    pred_corr(kk, 1) = Safe_Spearman(cands{kk, 2}, c_pred);
    pred_corr(kk, 2) = Safe_Spearman(cands{kk, 2}, difficulty);
    fprintf('  %-18s  vs c: %6.3f   vs difficulty (Newton steps after N_c=1): %6.3f\n', cands{kk, 1}, pred_corr(kk, 1), pred_corr(kk, 2));
end
fprintf('(If eta*||beta_dot|| correlates better than eta, the curvature limit lacks a path-length factor.)\n');

results = struct();
results.settings = struct('num_data_points', num_data_points, 'num_samples', S, 'fixed_sweep', fixed_sweep, ...
                          'num_max_newton_steps', num_max_newton_steps, 'newton_gradient_tol', newton_gradient_tol, ...
                          'newton_cost', newton_cost, 'eps_nl_both', eps_nl_both, 'eps_curv_sweep', eps_curv_sweep, ...
                          'pred_dts', pred_dts);
results.configs = configs;
results.fixed = fixed;  results.adapt = adapt;
results.best_fixed_index = jb;  results.oracle_cost = oracle_cost;  results.N_opt = N_opt;
results.pred = pred;  results.pred_corr = pred_corr;
results.eta = eta_s;  results.A = A_s;  results.difficulty = difficulty;
results_file = fullfile(ci_dir, 'ci-curvature-sweep-results.mat');
save(results_file, 'results');
fprintf('\nSaved results to:\n  %s\n', results_file);

Plot_CI_Curvature_Sweep;

%% Local helpers
function Print_Summary(s)
    fprintf('%-28s |  %4.2f / %2d   |   %4.2f / %2d     |   %5.2f / %5.2f / %5.2f |   %5.1f%%  | %.2f\n', ...
            s.label, s.mean_cont, max(s.cont_steps), s.mean_newton, max(s.newton_steps), ...
            s.mean_cost, s.p90_cost, s.max_cost, 100 * s.converged_frac, s.time_per_sample);
end

function rho = Safe_Spearman(x, y)
    ok = isfinite(x(:)) & isfinite(y(:));
    if nnz(ok) < 3 || numel(unique(x(ok))) < 2 || numel(unique(y(ok))) < 2
        rho = NaN;
    else
        rho = corr(x(ok), y(ok), 'Type', 'Spearman');
    end
end
