%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plots for Experiment 1 (Driver_CI_Mean_To_Sample), read from its saved
% results.  Bare-bones on purpose: plot/semilogy, labels, title, legend.

results_file = fullfile(fileparts(mfilename('fullpath')), 'ci-mean-to-sample-results.mat');
loaded = load(results_file, 'results');
res = loaded.results;

steps = res.settings.step_counts;
exp1 = res.exp1;
num_samples = res.settings.num_samples;

nom_err = arrayfun(@(e) median(e.nominal.ref_err, 'omitnan'), exp1);
mts_err = arrayfun(@(e) median(e.mts.ref_err, 'omitnan'), exp1);
nom_grad = arrayfun(@(e) median(e.nominal.grad_norm, 'omitnan'), exp1);
mts_grad = arrayfun(@(e) median(e.mts.grad_norm, 'omitnan'), exp1);
nom_time = arrayfun(@(e) e.nominal.time, exp1) / num_samples;
% The posterior mean is assumed known, so the mean solve is not charged.
mts_time = arrayfun(@(e) e.mts.time, exp1) / num_samples;

% Mean-to-sample + M extra Newton steps, one curve per M (absent in results
% saved before these variants were added; older files hold a single polish_M).
leg = {'nominal-to-sample', 'mean-to-sample'};
newton_err = [];  newton_grad = [];  newton_time = [];
if isfield(exp1, 'mts_newton')
    if isfield(res.settings, 'polish_Ms')
        polish_Ms = res.settings.polish_Ms;
        runs_of = @(e, kk) e.mts_newton{kk};
    else
        polish_Ms = res.settings.polish_M;
        runs_of = @(e, kk) e.mts_newton;
    end
    for kk = 1:length(polish_Ms)
        newton_err(:, kk) = arrayfun(@(e) median(runs_of(e, 3-kk).ref_err, 'omitnan'), exp1);
        newton_grad(:, kk) = arrayfun(@(e) median(runs_of(e, 3-kk).grad_norm, 'omitnan'), exp1);
        newton_time(:, kk) = arrayfun(@(e) runs_of(e, 3-kk).time, exp1) / num_samples;
        leg{end + 1} = sprintf('mean-to-sample + %d Newton', polish_Ms(3-kk)); %#ok<SAGROW>
    end
end
has_newton = ~isempty(newton_err);
newton_markers = {'^-', 'v-', 'd-', 'p-'};

% Reference: gradient norm of the low-fidelity objective at the low-fidelity
% optimizer (absent in results saved before it was added).
has_ref = isfield(res, 'lowfi_grad_norm');

%% 1. Accuracy vs steps
disp('Plot 1: accuracy vs steps');
figure;
semilogy(steps, nom_err, 'o-', steps, mts_err, 's-');
if has_newton
    hold on;
    for kk = 1:size(newton_err, 2), semilogy(steps, newton_err(:, kk), newton_markers{mod(kk - 1, 4) + 1}); end
    hold off;
end
xlabel('Continuation steps');
ylabel('Median relative error');
title('Accuracy vs steps');
legend(leg);
drawnow;

%% 2. Endpoint gradient vs steps
disp('Plot 2: endpoint gradient vs steps');
figure;
semilogy(steps, nom_grad, 'o-', steps, mts_grad, 's-');
if has_newton
    hold on;
    for kk = 1:size(newton_grad, 2), semilogy(steps, newton_grad(:, kk), newton_markers{mod(kk - 1, 4) + 1}); end
    hold off;
end
xlabel('Continuation steps');
ylabel('Median endpoint gradient norm');
if has_ref, hold on; yline(res.lowfi_grad_norm, 'k--'); hold off; end
title('Optimality vs steps');
if has_ref, legend([leg, {'low-fidelity optimizer'}]); else, legend(leg); end
drawnow;

%% 3. Accuracy vs time
disp('Plot 3: accuracy vs time');
figure;
semilogy(nom_time, nom_err, 'o-', mts_time, mts_err, 's-');
if has_newton
    hold on;
    for kk = 1:size(newton_err, 2), semilogy(newton_time(:, kk), newton_err(:, kk), newton_markers{mod(kk - 1, 4) + 1}); end
    hold off;
end
xlabel('Time per sample (s)');
ylabel('Median relative error');
title('Accuracy vs time');
legend(leg, 'Location', 'northeastoutside');
drawnow;

%% 4. Endpoint gradient vs time
disp('Plot 4: endpoint gradient vs time');
figure;
semilogy(nom_time, nom_grad, 'o-', mts_time, mts_grad, 's-');
if has_newton
    hold on;
    for kk = 1:size(newton_grad, 2), semilogy(newton_time(:, kk), newton_grad(:, kk), newton_markers{mod(kk - 1, 4) + 1}); end
    hold off;
end
if has_ref, hold on; yline(res.lowfi_grad_norm, 'k--', 'LineWidth', 2); hold off; end
xlabel('Time per sample (s)');
ylabel('Median endpoint gradient norm');
title('Optimality vs time');
if has_ref, legend([leg, {'low-fidelity optimizer'}], 'Location', 'northeastoutside'); else, legend(leg, 'Location', 'northeastoutside'); end
drawnow;

disp('Mean-to-sample plots done.');
