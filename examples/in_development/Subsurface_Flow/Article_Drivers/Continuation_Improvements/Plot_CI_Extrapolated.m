%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plots for Experiment 5 (Driver_CI_Extrapolated), read from its saved results:
% converged fraction and mean cost against the increment inflation lambda, one
% curve per method (fixed N_c, backtracking only, a priori heuristic).

results_file = fullfile(fileparts(mfilename('fullpath')), 'ci-extrapolated-results.mat');
loaded = load(results_file, 'results');
res = loaded.results;
S = res.settings;
nM = length(S.method_labels);
nL = length(S.lambdas);

conv = zeros(nM, nL);  cost = zeros(nM, nL);  steps = zeros(nM, nL);
for mm = 1:nM
    for ll = 1:nL
        s = res.summaries{mm, ll};
        conv(mm, ll) = 100 * s.converged_frac;
        cost(mm, ll) = s.mean_cost;
        steps(mm, ll) = s.mean_cont;
    end
end
nfix = length(S.fixed_steps);
heur = find(strcmp(S.method_labels, 'adaptive heuristic'));
bt_idx = find(strcmp(S.method_labels, 'backtrack only'));   % empty if it was not run
fixed_colors = lines(max(nfix, 1));
styles = repmat({'-o'}, 1, nM);
styles{heur} = 'r-d';
if ~isempty(bt_idx), styles{bt_idx} = 'k--s'; end

% Cost counted only where every sample converged (otherwise NaN: the method failed).
cost_ok = cost;
cost_ok(conv < 100 - 1.e-9) = NaN;
% Cheapest fixed N_c that fully converges, per lambda (the hindsight oracle).
oracle = min(cost_ok(1:nfix, :), [], 1);

figure; hold on;
for mm = 1:nM
    semilogx(S.lambdas, conv(mm, :), styles{mm});
end
set(gca, 'XScale', 'log');
xlabel('Increment scale \lambda'); ylabel('Samples converged (%)');
% legend(S.method_labels, 'Location', 'southwest'); 
grid on; drawnow;

figure; hold on;
for mm = 1:nM
    loglog(S.lambdas, cost_ok(mm, :), styles{mm});
end
loglog(S.lambdas, oracle, 'g-', 'LineWidth', 5);
xlabel('Increment inflation \lambda'); ylabel('Mean cost (fully converged runs only)');
title('Cost'); legend([S.method_labels, {'best fixed (hindsight)'}], 'Location', 'northwest'); grid on; drawnow;

figure; hold on;
for mm = 1:nM
    loglog(S.lambdas, steps(mm, :), styles{mm});
end
xlabel('Increment scale \lambda'); ylabel('Mean continuation steps');
legend(S.method_labels, 'Location', 'northeastoutside'); grid on; drawnow;

if false
    % Single summary figure: time per CONVERGED sample against lambda.
    %   time per converged sample = (batch time per sample) / (fraction converged)
    % A method that fails samples must pay for them, so its curve climbs (and ends
    % if it converges nowhere); the best method stays lowest.  Open markers = runs
    % that did not converge on every sample, labelled with the percentage.
    % Set cost_metric = 'work' to use steps + 0.5*Newton steps per converged sample
    % instead of wall-clock time (the batch time is parallel and noisy).
    cost_metric = 'time';

    metric = zeros(nM, nL);
    for mm = 1:nM
        for ll = 1:nL
            sm = res.summaries{mm, ll};
            if strcmp(cost_metric, 'work')
                per_sample = sm.mean_cost;
            else
                per_sample = sm.time_per_sample;
            end
            metric(mm, ll) = per_sample / max(sm.converged_frac, eps);
            if sm.converged_frac == 0
                metric(mm, ll) = NaN;
            end
        end
    end
    if strcmp(cost_metric, 'work')
        ylab = 'Cost per converged sample (steps + 0.5 Newton)';
    else
        ylab = 'Time per converged sample (s)';
    end

    h_line = gobjects(nM, 1);
    figure; hold on;
    for mm = 1:nM
        if mm == heur
            col = [0.84 0.15 0.16];  lw = 4;  ms = 12;           % highlighted
        elseif ~isempty(bt_idx) && mm == bt_idx
            col = [0.12 0.47 0.71];  lw = 2;  ms = 9;            % backtracking only
        else
            col = fixed_colors(mod(mm - 1, size(fixed_colors, 1)) + 1, :);  lw = 2;  ms = 8;
        end
        ls = '-';
        if mm <= nfix, ls = '--'; end
        h_line(mm) = plot(S.lambdas, metric(mm, :), ls, 'Color', col, 'LineWidth', lw);
        for ll = 1:nL
            if isnan(metric(mm, ll)), continue; end
            full = conv(mm, ll) >= 100 - 1.e-9;
            if full, mface = col; else, mface = 'w'; end
            plot(S.lambdas(ll), metric(mm, ll), 'o', 'Color', col, 'MarkerFaceColor', mface, ...
                'MarkerSize', ms, 'LineWidth', 2);
            if ~full
                text(S.lambdas(ll), metric(mm, ll), sprintf('  %.0f%%', conv(mm, ll)), 'Color', col, ...
                    'FontSize', 12, 'VerticalAlignment', 'bottom');
            end
        end
    end
    set(gca, 'XScale', 'log', 'YScale', 'log');
    xticks(S.lambdas);
    xlabel('Increment inflation \lambda (harder to the right)');
    ylabel(ylab);
    title('Cost to obtain a converged sample (open marker: not all samples converged)');
    legend(h_line, S.method_labels, 'Location', 'northwest'); grid on; drawnow;
end
disp('Extrapolated-increment plots done.');
