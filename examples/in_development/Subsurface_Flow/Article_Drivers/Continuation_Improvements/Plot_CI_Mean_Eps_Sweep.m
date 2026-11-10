%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plots for Experiment 6 (Driver_CI_Mean_Eps_Sweep), read from its saved results.
% Figure 1: endpoint gradient norm of the mean continuation against continuation
% steps (+ backtracks) for the eps_nl sweep, with 0, 1 and 2 Newton polishing
% steps.  Points that share a step count lie on one vertical line; the eps_nl
% value(s) that produce it are written along the bottom of that line, in a band
% below the data so the labels never touch the curves.
% Figures 2 and 3 (eps_curv sweep; accuracy vs work with fixed N_c) are
% commented out below.

results_file = fullfile(fileparts(mfilename('fullpath')), 'ci-mean-eps-sweep-results.mat');
loaded = load(results_file, 'results');
res = loaded.results;
S = res.settings;

%% Figure 1
num_polish = S.num_polish;
nR = length(res.nl_runs);
w_nl = cellfun(@(x) x.steps + x.backtracks, res.nl_runs);
g_nl = zeros(num_polish + 1, nR);       % row m + 1: gradient norm after m Newton steps
for ii = 1:nR
    g_nl(:, ii) = res.nl_runs{ii}.newton_grad_norms(1:num_polish + 1).';
end

figure; hold on;
styles = {'-o', '-s', '-^', '-d'};
leg = cell(1, num_polish + 1);
for m = 0:num_polish
    semilogy(w_nl, g_nl(m + 1, :), styles{m + 1});   % one point per eps_nl, in sweep order
    leg{m + 1} = sprintf('%d Newton step%s', m, repmat('s', 1, m ~= 1));
end
yline(res.lowfi_grad_norm, 'k--', 'LineWidth', 4);
leg{end + 1} = 'low-fidelity optimizer';
set(gca, 'XScale', 'log', 'YScale', 'log');

% Reserve a label band below the data (a fixed share of the log-height).
y_all = [g_nl(:); res.lowfi_grad_norm];
y_all = y_all(y_all > 0 & isfinite(y_all));
lo = log10(min(y_all));  hi = log10(max(y_all));
span = 0.5*max(hi - lo, 1);
y_bottom = 10^(lo - 0.9 * span);
ylim([y_bottom, 10^(hi + 0.1 * span)]);
x_unique = unique(w_nl);
xlim([min(x_unique) / 1.3, max(x_unique) * 1.3]);

for x = x_unique(:).'
    eps_here = S.eps_nl_values(w_nl == x);
    parts = arrayfun(@(e) eps_text(e), eps_here, 'UniformOutput', false);
    lab = ['$\epsilon_{nl}=', strjoin(parts, ', '), '$'];
    plot([x x], [y_bottom, 10^(hi + 0.1 * span)], ':', 'Color', [0.6 0.6 0.6], 'HandleVisibility', 'off');
    text(x, y_bottom * 1.15, lab, 'Interpreter', 'latex', 'FontSize', 12, 'Rotation', 90, ...
         'HorizontalAlignment', 'left', 'VerticalAlignment', 'middle', 'BackgroundColor', 'w');
end

xlabel('Continuation steps'); ylabel('Endpoint gradient norm');
% title('Mean: nonlinearity limit only (\epsilon_{curv} = Inf)');
legend(leg, 'Location', 'northeastoutside'); grid on; drawnow;

%% Figures 2 and 3 (disabled for now)
% grad_of = @(runs) cellfun(@(x) x.grad_norm, runs);
% work_of = @(runs) cellfun(@(x) x.steps + x.backtracks, runs);
% g_cv = grad_of(res.curv_runs);    w_cv = work_of(res.curv_runs);
% g_fx = grad_of(res.fixed_runs);   w_fx = S.fixed_steps;
%
% figure; hold on;
% x_cv = S.eps_curv_values;  x_cv(isinf(x_cv)) = 10 * max(x_cv(isfinite(x_cv)));
% loglog(x_cv, g_cv, '-o');
% yline(res.lowfi_grad_norm, 'k--');
% set(gca, 'XScale', 'log', 'YScale', 'log', 'XDir', 'reverse');
% xlabel('\epsilon_{curv} (Inf = limit off, drawn at 10x the largest value)'); ylabel('||grad|| after continuation (no Newton)');
% title('Mean: curvature limit only (\epsilon_{nl} = Inf)');
% legend({'adaptive', 'low-fidelity optimizer'}, 'Location', 'best'); grid on; drawnow;
%
% figure; hold on;
% loglog(w_nl, g_nl(1, :), '-o');       % eps_nl sweep, no Newton
% loglog(w_cv, g_cv, '-s');
% loglog(w_fx, g_fx, 'k-^');
% yline(res.lowfi_grad_norm, 'k--');
% set(gca, 'XScale', 'log', 'YScale', 'log');
% xlabel('Continuation steps (+ backtracks)'); ylabel('||grad|| after continuation (no Newton)');
% title('Mean: accuracy vs work');
% legend({'eps_{nl} sweep', 'eps_{curv} sweep', 'fixed N_c', 'low-fidelity optimizer'}, 'Location', 'best'); grid on; drawnow;

disp('Mean eps-sweep plots done.');

function t = eps_text(e)
    if isinf(e)
        t = '\infty';
    else
        t = sprintf('%g', e);
    end
end
