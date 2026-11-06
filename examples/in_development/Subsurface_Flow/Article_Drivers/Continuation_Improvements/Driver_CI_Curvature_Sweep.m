%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Experiment 2c: is the curvature limit (eps_curv) of the a priori time-step
% heuristic needed, and how should it be set?
%
% Part 1 (ablation / frontier).  Mean-to-sample continuation is run with
%   - fixed N_c steps (the reference frontier),
%   - 'both'      : eps_nl = eps_nl_both,  eps_curv swept (Inf = nonlinearity limit only),
%   - 'curv only' : eps_nl = Inf,          eps_curv swept (curvature limit only),
% each with M = 0 and M = 2 Newton polishing steps (M = 0 shows what
% polishing would otherwise hide).  Work = steps + backtracks + M.  Each run
% is scored by the endpoint gradient norm ||g_beta|| (median, 90th percentile,
% max over samples), the fraction of samples reaching grad_target, and the
% ratio to the fixed-step frontier interpolated at the same mean work
% (ratio < 1: better than fixed steps at equal work).
%
% Part 2 (predictor-error test).  For each sample, one Euler predictor step
% from the mean start at dt = 0.5 and 0.25.  Its gradient defect
% ||g(pred) - g(start)|| should scale like dt^2; the coefficient c is compared
% (Spearman) with eta, eta * ||beta_dot||, A, A * ||beta_dot|| and ||beta_dot||.
% If eta * ||beta_dot|| predicts c better than eta, the curvature limit lacks
% a path-length factor.

ci_num_data_points = 1;                % this experiment uses 1 data point (other drivers use 3)
CI_Setup;
CI_Sample_Reference;
clear ci_num_data_points;

fixed_sweep = [1, 2, 3, 4, 5];
polish_Ms = [0, 2];                    % Newton polishing steps at t = 1 (counted as work)
eps_nl_both = 0.01;                    % eps_nl for the 'both' family
eps_curv_sweep = [1, 3, 10, 30, Inf];  % Inf = curvature limit off
grad_target = lowfi_grad_norm;         % "converged" level for the meet-target fractions
pred_dts = [0.5, 0.25];

% Adaptive configurations.
configs = struct('family', {}, 'eps_nl', {}, 'eps_curv', {});
for ec = eps_curv_sweep
    configs(end + 1) = struct('family', 'both', 'eps_nl', eps_nl_both, 'eps_curv', ec); %#ok<SAGROW>
end
for ec = eps_curv_sweep(isfinite(eps_curv_sweep))
    configs(end + 1) = struct('family', 'curv only', 'eps_nl', Inf, 'eps_curv', ec); %#ok<SAGROW>
end
num_cfg = length(configs);
num_M = length(polish_Ms);
num_N = length(fixed_sweep);
S = num_samples;

fprintf('\n============================================================\n');
fprintf('Experiment 2c: curvature-limit sweep (%d data points, %d samples)\n', num_data_points, S);
fprintf('============================================================\n');
fprintf('grad_target (low-fidelity optimizer ||grad||) = %.3e\n', grad_target);
fprintf('Work = steps + backtracks + Newton steps M.  "ratio" = adaptive stat / fixed-step\n');
fprintf('frontier stat interpolated at the same mean work (<1 means adaptive is better).\n');

fixed = cell(num_M, 1);
adapt = cell(num_M, 1);
for ip = 1:num_M
    M = polish_Ms(ip);
    fprintf('\n------------------------------------------------------------\n');
    fprintf('Newton polishing steps M = %d\n', M);
    fprintf('------------------------------------------------------------\n');

    fprintf('Fixed steps (mean-to-sample):\n');
    fprintf('  N_c | work | median grad | p90 grad    | max grad    | meet target | s/sample\n');
    fixed{ip} = [];
    for jj = 1:num_N
        run_out = Run_Samples(sen_op, md_hessian_analysis, z_opt, r, 'fixed_mts', fixed_sweep(jj), beta_bar, [], [], sample_indices, M);
        sm = Summarize_CI_Run(run_out, grad_target);
        sm.N_c = fixed_sweep(jj);
        fixed{ip} = [fixed{ip}, sm];
        fprintf('  %3d | %4.1f | %.4e  | %.4e  | %.4e  |   %5.1f%%    | %.2f\n', sm.N_c, sm.work_mean, ...
                sm.median, sm.p90, sm.max, 100 * sm.meet_frac, sm.time_per_sample);
    end
    fixed_work = [fixed{ip}.work_mean];

    fprintf('Adaptive (eps_curv = Inf means curvature limit off; eps_nl = Inf means nonlinearity limit off):\n');
    fprintf('  family    | eps_nl | eps_curv | work mean/max | steps | bt | curv binds | median (ratio)      | p90 (ratio)         | max (ratio)         | meet   | s/sample\n');
    adapt{ip} = [];
    for cc = 1:num_cfg
        cfg = configs(cc);
        bt_options = sample_bt_options;
        bt_options.epsilon_curv = cfg.eps_curv;
        run_out = Run_Samples(sen_op, md_hessian_analysis, z_opt, r, 'adaptive_mts', [], beta_bar, cfg.eps_nl, bt_options, sample_indices, M);
        sm = Summarize_CI_Run(run_out, grad_target);
        sm.family = cfg.family;  sm.eps_nl = cfg.eps_nl;  sm.eps_curv = cfg.eps_curv;
        sm.steps = run_out.steps;  sm.dt_prior = run_out.dt_prior;
        sm.A = run_out.A;  sm.eta = run_out.eta;  sm.dt_nl = run_out.dt_nl;  sm.dt_curv = run_out.dt_curv;
        sm.curv_binds = mean(sm.dt_curv < sm.dt_nl, 'omitnan');   % fraction of samples limited by curvature
        sm.ratio_median = Frontier_Ratio(fixed_work, [fixed{ip}.median], sm.work_mean, sm.median);
        sm.ratio_p90 = Frontier_Ratio(fixed_work, [fixed{ip}.p90], sm.work_mean, sm.p90);
        sm.ratio_max = Frontier_Ratio(fixed_work, [fixed{ip}.max], sm.work_mean, sm.max);
        adapt{ip} = [adapt{ip}, sm];
        fprintf('  %-9s | %6.3g | %8.3g | %5.1f / %2d    | %5.2f | %2d |   %5.1f%%   | %.3e (%s) | %.3e (%s) | %.3e (%s) | %5.1f%% | %.2f\n', ...
                sm.family, sm.eps_nl, sm.eps_curv, sm.work_mean, sm.work_max, sm.steps_mean, sm.backtracks, 100 * sm.curv_binds, ...
                sm.median, Ratio_String(sm.ratio_median), sm.p90, Ratio_String(sm.ratio_p90), sm.max, Ratio_String(sm.ratio_max), ...
                100 * sm.meet_frac, sm.time_per_sample);
    end
end

%% Smallest fixed N_c that reaches grad_target, per sample
N_star = nan(S, num_M);
for ip = 1:num_M
    for ii = 1:S
        for jj = 1:num_N
            if fixed{ip}(jj).grad_norm(ii) <= grad_target
                N_star(ii, ip) = fixed_sweep(jj);
                break;
            end
        end
    end
end

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

eta_s = adapt{1}(1).eta;               % a priori quantities do not depend on eps or M
A_s = adapt{1}(1).A;
one_step_grad = fixed{1}(1).grad_norm;  % unpolished one-step gradient: difficulty proxy
c_pred = pred.c(:, end);

fprintf('\nPer-sample table:\n');
fprintf(' sample |    A       |   eta      | ||beta_dot|| | c (pred)   | slope | 1-step grad |');
for ip = 1:num_M
    fprintf(' N* (M=%d) |', polish_Ms(ip));
end
fprintf('\n');
for ii = 1:S
    fprintf(' %5d  | %.3e | %.3e | %.3e   | %.3e | %5.2f | %.3e   |', ii, A_s(ii), eta_s(ii), pred.speed(ii), c_pred(ii), pred.slope(ii), one_step_grad(ii));
    for ip = 1:num_M
        if isnan(N_star(ii, ip))
            fprintf('   >%d     |', max(fixed_sweep));
        else
            fprintf('   %3d     |', N_star(ii, ip));
        end
    end
    fprintf('\n');
end
fprintf('(N* = smallest fixed N_c reaching grad_target; slope should be ~2 for a clean quadratic predictor defect.)\n');

fprintf('\nSpearman correlation with the predictor defect coefficient c (and with the 1-step gradient):\n');
cands = {'eta', eta_s; 'eta*||beta_dot||', eta_s .* pred.speed; 'A', A_s; 'A*||beta_dot||', A_s .* pred.speed; '||beta_dot||', pred.speed};
pred_corr = zeros(size(cands, 1), 2);
for kk = 1:size(cands, 1)
    pred_corr(kk, 1) = Safe_Spearman(cands{kk, 2}, c_pred);
    pred_corr(kk, 2) = Safe_Spearman(cands{kk, 2}, one_step_grad);
    fprintf('  %-18s  vs c: %6.3f   vs 1-step grad: %6.3f\n', cands{kk, 1}, pred_corr(kk, 1), pred_corr(kk, 2));
end
fprintf('(If eta*||beta_dot|| correlates better than eta, the curvature limit lacks a path-length factor.)\n');

results = struct();
results.settings = struct('num_data_points', num_data_points, 'num_samples', S, 'fixed_sweep', fixed_sweep, ...
                          'polish_Ms', polish_Ms, 'eps_nl_both', eps_nl_both, 'eps_curv_sweep', eps_curv_sweep, ...
                          'grad_target', grad_target, 'pred_dts', pred_dts);
results.configs = configs;
results.fixed = fixed;  results.adapt = adapt;  results.N_star = N_star;
results.pred = pred;  results.pred_corr = pred_corr;
results.eta = eta_s;  results.A = A_s;  results.one_step_grad = one_step_grad;
results_file = fullfile(ci_dir, 'ci-curvature-sweep-results.mat');
save(results_file, 'results');
fprintf('\nSaved results to:\n  %s\n', results_file);

Plot_CI_Curvature_Sweep;

%% Local helpers
function s = Summarize_CI_Run(run_out, grad_target)
    g = run_out.grad_norm(:);
    ok = isfinite(g);
    s.grad_norm = g;
    s.work_mean = mean(run_out.work, 'omitnan');
    s.work_max = max(run_out.work, [], 'omitnan');
    s.steps_mean = mean(run_out.steps, 'omitnan');
    s.backtracks = sum(run_out.backtracks);
    s.median = median(g, 'omitnan');
    s.max = max(g, [], 'omitnan');
    if any(ok)
        s.p90 = prctile(g(ok), 90);
        s.meet_frac = mean(g(ok) <= grad_target);
    else
        s.p90 = NaN;
        s.meet_frac = NaN;
    end
    s.failures = sum(run_out.failures);
    s.time_per_sample = run_out.time / numel(g);
end

function ratio = Frontier_Ratio(fixed_work, fixed_stat, work, stat)
    % adaptive stat / fixed-step stat interpolated (log-linearly) at the same work
    ratio = NaN;
    if work < min(fixed_work) || work > max(fixed_work) || any(~isfinite(fixed_stat)) || any(fixed_stat <= 0)
        return;
    end
    ratio = stat / 10^interp1(fixed_work, log10(fixed_stat), work);
end

function str = Ratio_String(ratio)
    if isnan(ratio)
        str = '  n/a';
    else
        str = sprintf('%5.2f', ratio);
    end
end

function rho = Safe_Spearman(x, y)
    ok = isfinite(x(:)) & isfinite(y(:));
    if nnz(ok) < 3 || numel(unique(x(ok))) < 2 || numel(unique(y(ok))) < 2
        rho = NaN;
    else
        rho = corr(x(ok), y(ok), 'Type', 'Spearman');
    end
end
