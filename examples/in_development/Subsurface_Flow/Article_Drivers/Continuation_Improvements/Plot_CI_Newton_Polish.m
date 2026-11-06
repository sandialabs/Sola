% Plots for Experiment 4 (bare-bones): median and max endpoint gradient norm
% against continuation steps (one curve per number of extra Newton steps) and
% against time per sample.
loaded = load(fullfile(fileparts(mfilename('fullpath')), 'ci-newton-polish-results.mat'));
res = loaded.results;
S = res.settings;
method = S.methods{1};
runs = res.exp4.(method);
nN = length(S.step_counts);
nM = length(S.newton_counts);

med_grad = zeros(nN, nM);  max_grad = zeros(nN, nM);  t_per = zeros(nN, nM);
for ii = 1:nN
    for jj = 1:nM
        med_grad(ii, jj) = median(runs{ii, jj}.grad_norm, 'omitnan');
        max_grad(ii, jj) = max(runs{ii, jj}.grad_norm);
        t_per(ii, jj) = runs{ii, jj}.time / S.num_samples;
    end
end
leg = arrayfun(@(m) sprintf('%d Newton Steps', m), S.newton_counts, 'UniformOutput', false);

% Reference: gradient norm of the low-fidelity objective at the low-fidelity
% optimizer; only present in results saved by newer versions of the driver.
has_ref = isfield(res, 'lowfi_grad_norm');
ref_name = 'Low-fidelity optimizer';

figure(1);
semilogy(S.step_counts, med_grad, '-o');
if has_ref, hold on; yline(res.lowfi_grad_norm, 'k--', ref_name); hold off; end
xlabel('Continuation steps N_c'); ylabel('Median gradient norm');
title('Median endpoint gradient vs continuation steps'); legend(leg); grid on; drawnow;

figure(2);
semilogy(S.step_counts, max_grad, '-o');
if has_ref, hold on; yline(res.lowfi_grad_norm, 'k--', ref_name); hold off; end
xlabel('Continuation steps N_c'); ylabel('Max gradient norm');
title('Max endpoint gradient vs continuation steps'); legend(leg); grid on; drawnow;

figure(3);
loglog(t_per, med_grad, '-o');
if has_ref, hold on; yline(res.lowfi_grad_norm, 'k--', ref_name, 'LabelHorizontalAlignment', 'left', 'LabelVerticalAlignment', 'bottom', 'LineWidth', 2); hold off; end
xlabel('Time per sample (s)'); ylabel('Median gradient norm');
legend(leg); grid on; drawnow;
