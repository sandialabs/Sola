%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Experiment 2d: how many Newton polishing steps, and is it better to polish
% than to take more continuation steps?
%
% Mean-to-sample continuation is run with a few step counts (fixed N_c and two
% adaptive settings) and M_max Newton steps at t = 1.  The Newton iterates are
% prefixes of one another, so a single run gives the endpoint gradient norm
% ||g_beta|| after every M = 0..M_max.  This lets us compare
%   - fixed M (every sample gets M Newton steps), and
%   - a stopping rule: Newton until ||g_beta|| <= tol (per sample), capped at
%     M_max.  The gradient is already computed by each Newton step, so the
%     rule is free; its cost is the per-sample number of steps used.
% Cost = continuation steps (+ backtracks) + newton_cost * Newton steps.
% newton_cost = 0.5 matches the earlier timings (a Newton step takes about half
% as long as a continuation step: one gradient + one solve vs. two gradients,
% B and two solves).  The Pareto-optimal (cost, gradient) settings are listed.

ci_num_data_points = 1;                % this experiment uses 1 data point (other drivers use 3)
CI_Setup;
CI_Sample_Reference;
clear ci_num_data_points;

fixed_Ncs = [1, 2, 3];                 % fixed continuation step counts
adaptive_eps_nl = 0.01;
adaptive_eps_curvs = [3, 10];          % adaptive settings (eps_nl, eps_curv)
M_max = 6;                             % Newton steps run; all M = 0..M_max are scored
newton_cost = 0.5;                     % cost of a Newton step relative to a continuation step
grad_target = lowfi_grad_norm;         % "converged" level
stop_tols = grad_target * [1, 10, 100];  % Newton stopping tolerances

S = num_samples;
num_fixed = length(fixed_Ncs);
num_adapt = length(adaptive_eps_curvs);
num_cfg = num_fixed + num_adapt;

fprintf('\n============================================================\n');
fprintf('Experiment 2d: Newton polishing sweep (%d data points, %d samples)\n', num_data_points, S);
fprintf('============================================================\n');
fprintf('grad_target (low-fidelity optimizer ||grad||) = %.3e\n', grad_target);
fprintf('cost = continuation steps (+ backtracks) + %.2f * Newton steps\n', newton_cost);

%% Run
cfgs = struct('name', {}, 'steps', {}, 'gn', {}, 'time_per_sample', {});
for jj = 1:num_cfg
    if jj <= num_fixed
        name = sprintf('fixed N_c = %d', fixed_Ncs(jj));
        run_out = Run_Samples(sen_op, md_hessian_analysis, z_opt, r, 'fixed_mts', fixed_Ncs(jj), beta_bar, [], [], sample_indices, M_max);
        steps = fixed_Ncs(jj) * ones(S, 1);
    else
        ec = adaptive_eps_curvs(jj - num_fixed);
        name = sprintf('adaptive eps_nl = %g, eps_curv = %g', adaptive_eps_nl, ec);
        bt_options = sample_bt_options;
        bt_options.epsilon_curv = ec;
        run_out = Run_Samples(sen_op, md_hessian_analysis, z_opt, r, 'adaptive_mts', [], beta_bar, adaptive_eps_nl, bt_options, sample_indices, M_max);
        steps = run_out.steps + run_out.backtracks;
    end
    cfgs(jj).name = name;
    cfgs(jj).steps = steps;
    cfgs(jj).gn = run_out.newton_grad_norms;
    cfgs(jj).time_per_sample = run_out.time / S;
    fprintf('ran %-40s | mean steps %.2f | %.2f s/sample (incl. %d Newton steps)\n', name, mean(steps, 'omitnan'), run_out.time / S, M_max);
end

%% Score every (configuration, rule)
points = struct('label', {}, 'cfg', {}, 'kind', {}, 'cost', {}, 'median', {}, 'p90', {}, 'max', {}, 'meet', {});
for jj = 1:num_cfg
    gn = cfgs(jj).gn;
    steps = cfgs(jj).steps;
    mean_steps = mean(steps, 'omitnan');

    fprintf('\n------------------------------------------------------------\n');
    fprintf('%s   (mean continuation steps %.2f)\n', cfgs(jj).name, mean_steps);
    fprintf('------------------------------------------------------------\n');

    % Per-step reduction factor of the gradient (quadratic convergence would accelerate).
    red = nan(1, M_max);
    for m = 1:M_max
        red(m) = median(gn(:, m + 1) ./ gn(:, m), 'omitnan');
    end
    fprintf('Median gradient reduction factor g(M+1)/g(M) for M = 0..%d:\n  %s\n', M_max - 1, sprintf('%.2e  ', red));

    fprintf('  rule              | Newton M mean/max | cost  | median      | p90         | max         | meet target\n');
    fixedM = struct('m', {}, 'cost', {}, 'median', {}, 'p90', {}, 'max', {}, 'meet', {});
    for m = 0:M_max
        g = gn(:, m + 1);
        st = Score_Gradients(g, grad_target);
        cost = mean_steps + newton_cost * m;
        fixedM(end + 1) = struct('m', m, 'cost', cost, 'median', st.median, 'p90', st.p90, 'max', st.max, 'meet', st.meet); %#ok<SAGROW>
        points(end + 1) = struct('label', sprintf('%s, M=%d', cfgs(jj).name, m), 'cfg', jj, 'kind', 'fixed M', ...
                                 'cost', cost, 'median', st.median, 'p90', st.p90, 'max', st.max, 'meet', st.meet); %#ok<SAGROW>
        fprintf('  fixed M = %d       |   %d / %d           | %5.2f | %.4e  | %.4e  | %.4e  |   %5.1f%%\n', ...
                m, m, m, cost, st.median, st.p90, st.max, 100 * st.meet);
    end

    stopM = struct('tol', {}, 'cost', {}, 'mean_M', {}, 'max_M', {}, 'median', {}, 'p90', {}, 'max', {}, 'meet', {}, 'reached', {});
    for tol = stop_tols
        m_used = nan(S, 1);
        g_final = nan(S, 1);
        reached = false(S, 1);
        for ii = 1:S
            row = gn(ii, :);
            if all(isnan(row))
                continue;
            end
            j = find(row <= tol, 1);
            if isempty(j)
                j = M_max + 1;
            else
                reached(ii) = true;
            end
            m_used(ii) = j - 1;
            g_final(ii) = row(j);
        end
        st = Score_Gradients(g_final, grad_target);
        cost = mean(steps + newton_cost * m_used, 'omitnan');
        stopM(end + 1) = struct('tol', tol, 'cost', cost, 'mean_M', mean(m_used, 'omitnan'), 'max_M', max(m_used), ...
                                'median', st.median, 'p90', st.p90, 'max', st.max, 'meet', st.meet, 'reached', mean(reached)); %#ok<SAGROW>
        points(end + 1) = struct('label', sprintf('%s, stop tol=%.2e', cfgs(jj).name, tol), 'cfg', jj, 'kind', 'stop', ...
                                 'cost', cost, 'median', st.median, 'p90', st.p90, 'max', st.max, 'meet', st.meet); %#ok<SAGROW>
        fprintf('  stop tol %.2e |   %.2f / %d         | %5.2f | %.4e  | %.4e  | %.4e  |   %5.1f%%  (tol reached: %.0f%%)\n', ...
                tol, mean(m_used, 'omitnan'), max(m_used), cost, st.median, st.p90, st.max, 100 * st.meet, 100 * mean(reached));
    end
    cfgs(jj).mean_steps = mean_steps;
    cfgs(jj).fixedM = fixedM;
    cfgs(jj).stopM = stopM;

    bad = find(gn(:, end) > grad_target | isnan(gn(:, end)));
    if isempty(bad)
        fprintf('All samples reach grad_target by M = %d.\n', M_max);
    else
        fprintf('Samples above grad_target at M = %d: %s\n', M_max, mat2str(bad'));
        fprintf('  their ||g_beta|| at M = %d: %s\n', M_max, sprintf('%.2e ', gn(bad, end)));
    end
end

%% Pareto-optimal settings
for stat = {'p90', 'median'}
    st_name = stat{1};
    vals = [points.(st_name)];
    costs = [points.cost];
    keep = Pareto_Indices(costs, vals);
    fprintf('\n============================================================\n');
    fprintf('Pareto-optimal settings in (cost, %s gradient), cheapest first:\n', st_name);
    fprintf('============================================================\n');
    fprintf('  cost  | %-10s | meet   | setting\n', st_name);
    for k = keep
        fprintf('  %5.2f | %.4e | %5.1f%% | %s\n', points(k).cost, vals(k), 100 * points(k).meet, points(k).label);
    end
end

results = struct();
results.settings = struct('num_data_points', num_data_points, 'num_samples', S, 'fixed_Ncs', fixed_Ncs, ...
                          'adaptive_eps_nl', adaptive_eps_nl, 'adaptive_eps_curvs', adaptive_eps_curvs, ...
                          'M_max', M_max, 'newton_cost', newton_cost, 'grad_target', grad_target, 'stop_tols', stop_tols);
results.cfgs = cfgs;  results.points = points;
results_file = fullfile(ci_dir, 'ci-polish-sweep-results.mat');
save(results_file, 'results');
fprintf('\nSaved results to:\n  %s\n', results_file);

Plot_CI_Polish_Sweep;

%% Local helpers
function st = Score_Gradients(g, grad_target)
    ok = isfinite(g);
    st.median = median(g, 'omitnan');
    st.max = max(g, [], 'omitnan');
    if any(ok)
        st.p90 = prctile(g(ok), 90);
        st.meet = mean(g(ok) <= grad_target);
    else
        st.p90 = NaN;
        st.meet = NaN;
    end
end

function keep = Pareto_Indices(cost, stat)
    % indices whose stat is strictly below that of every cheaper-or-equal-cost point
    [~, order] = sort(cost);
    best = inf;
    keep = [];
    for k = order(:)'
        if isfinite(stat(k)) && stat(k) < best
            keep(end + 1) = k; %#ok<AGROW>
            best = stat(k);
        end
    end
end
