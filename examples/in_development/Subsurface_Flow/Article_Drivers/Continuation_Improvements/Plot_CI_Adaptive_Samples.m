%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plots for Experiment 2b (Driver_CI_Adaptive_Samples), read from its saved
% results.  Bare-bones on purpose: plot/semilogy, labels, title, legend.

results_file = fullfile(fileparts(mfilename('fullpath')), 'ci-adaptive-samples-results.mat');
loaded = load(results_file, 'results');
res = loaded.results;

num_samples = res.settings.num_samples;
polish_M = res.settings.polish_M;      % Newton polishing steps at t = 1 (counted as work)
work_label = 'Work per sample (steps + backtracks + Newton steps)';
eps_list = res.settings.sample_epsilon_nl;
adapt = res.sample_adapt{res.settings.show_idx};    % single run shown in per-sample plots
adapt_label = sprintf('adaptive (\\epsilon_{nl} = %g)', eps_list(res.settings.show_idx));

%% 5. Adaptive steps, samples
disp('Plot 5: adaptive steps, samples');
fixed_err = cellfun(@(s) median(s.grad_norm, 'omitnan'), res.sample_fixed);
figure;
adapt_work_mean = cellfun(@(a) mean(a.work, 'omitnan'), res.sample_adapt);
adapt_err = cellfun(@(a) median(a.grad_norm, 'omitnan'), res.sample_adapt);
semilogy(res.settings.sample_fixed_sweep + res.settings.polish_M, fixed_err, 'o-', adapt_work_mean, adapt_err, 's-');
xlabel(work_label);
ylabel('Median endpoint gradient norm');
title(sprintf('Samples: fixed vs adaptive (%d Newton steps)', polish_M));
legend('fixed', 'adaptive');
drawnow;

%% 6. Adaptive steps, per-sample work
disp('Plot 6: per-sample work');
figure;
plot(1:num_samples, adapt.steps, 'o-', 1:num_samples, adapt.backtracks, 's-');
xlabel('Sample index');
ylabel('Count');
title(['Adaptive work per sample, ' adapt_label]);
legend('accepted steps', 'backtracks');
drawnow;

%% 7. Tail of the per-sample gradient vs mean work
% An allocator should help the worst samples, not the median one.  Failed
% adaptive samples have NaN gradients and are omitted from these statistics.
disp('Plot 7: gradient tail vs mean work');
fixed_work = res.settings.sample_fixed_sweep + res.settings.polish_M;
fixed_p90 = cellfun(@(s) prctile(s.grad_norm(isfinite(s.grad_norm)), 90), res.sample_fixed);
fixed_max = cellfun(@(s) max(s.grad_norm, [], 'omitnan'), res.sample_fixed);
adapt_work_mean = cellfun(@(a) mean(a.work, 'omitnan'), res.sample_adapt);
adapt_p90 = cellfun(@(a) prctile(a.grad_norm(isfinite(a.grad_norm)), 90), res.sample_adapt);
adapt_max = cellfun(@(a) max(a.grad_norm, [], 'omitnan'), res.sample_adapt);
figure;
semilogy(fixed_work, fixed_p90, 'o-', fixed_work, fixed_max, 'o--', ...
         adapt_work_mean, adapt_p90, 's-', adapt_work_mean, adapt_max, 's--');
xlabel(work_label);
ylabel('Endpoint gradient norm');
title('Worst-case samples: fixed vs adaptive');
legend('fixed (90th pct)', 'fixed (max)', 'adaptive (90th pct)', 'adaptive (max)');
drawnow;

%% 8. Per-sample endpoint gradient: fixed steps vs adaptive
% Samples are sorted by the number of steps the adaptive run chose.  If the
% heuristic assigns smaller steps to samples that need them, the fixed-step
% gradient norms should be largest toward the right (many adaptive steps).
disp('Plot 8: per-sample gradient, fixed vs adaptive');
compare_N = res.settings.sample_fixed_sweep;   % all fixed step counts that were run
adapt_work = adapt.work(:);
[~, order] = sort(adapt_work);
figure;
hold on;
legend_labels = cell(1, length(compare_N) + 1);
for kk = 1:length(compare_N)
    j = find(res.settings.sample_fixed_sweep == compare_N(kk), 1);
    semilogy(1:num_samples, res.sample_fixed{j}.grad_norm(order), 'o-');
    legend_labels{kk} = sprintf('fixed N_c = %d', compare_N(kk));
end
semilogy(1:num_samples, adapt.grad_norm(order), 'k s-');
legend_labels{end} = adapt_label;
set(gca, 'YScale', 'log');
xlabel('Sample (sorted by adaptive work)');
ylabel('Endpoint gradient norm');
title('Per-sample gradient: fixed vs adaptive');
legend(legend_labels);
drawnow;

%% 9. Does the heuristic's time step track sample difficulty?
% Difficulty proxy: gradient norm after a single fixed step.  One series per
% epsilon_nl; with the curvature limit switched off, dt should fall as the
% sample gets harder.
disp('Plot 9: heuristic time step vs sample difficulty');
j1 = find(res.settings.sample_fixed_sweep == 1, 1);
difficulty = res.sample_fixed{j1}.grad_norm(:);
figure;
hold on;
for kk = 1:length(res.sample_adapt)
    dt = res.sample_adapt{kk}.dt_prior(:);
    ok = isfinite(difficulty) & isfinite(dt);
    fprintf('eps_nl = %g: Spearman correlation (difficulty vs dt_prior) = %.3f\n', ...
            eps_list(kk), corr(difficulty(ok), dt(ok), 'Type', 'Spearman'));
    loglog(difficulty, dt, 'o');
end
set(gca, 'XScale', 'log', 'YScale', 'log');
xlabel('Gradient norm after one fixed step + Newton polishing');
ylabel('Heuristic time step dt');
title('A priori time step vs sample difficulty');
legend(arrayfun(@(e) sprintf('\\epsilon_{nl} = %g', e), eps_list, 'UniformOutput', false));
drawnow;

disp('Adaptive-samples plots done.');
