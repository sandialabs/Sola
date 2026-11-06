%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plots for Experiment 4b (Driver_CI_Mean_Newton_Polish), read from its saved
% results: posterior-mean error and endpoint gradient norm against continuation
% steps and time (one curve per number of extra Newton steps).

results_file = fullfile(fileparts(mfilename('fullpath')), 'ci-mean-newton-polish-results.mat');
loaded = load(results_file, 'results');
res = loaded.results;
S = res.settings;
nN = length(S.step_counts);
nM = length(S.newton_counts);

err = zeros(nN, nM);  grad = zeros(nN, nM);  t_run = zeros(nN, nM);
for ii = 1:nN
    for jj = 1:nM
        err(ii, jj) = res.runs{ii, jj}.ref_err;
        grad(ii, jj) = res.runs{ii, jj}.grad_norm;
        t_run(ii, jj) = res.runs{ii, jj}.time;
    end
end
leg = arrayfun(@(m) sprintf('M = %d Newton', m), S.newton_counts, 'UniformOutput', false);

figure;
semilogy(S.step_counts, err, '-o');
xlabel('Continuation steps N_c'); ylabel('Relative error');
title('Posterior mean: error vs steps'); legend(leg); grid on; drawnow;

figure;
semilogy(S.step_counts, grad, '-o');
yline(res.lowfi_grad_norm, 'k--');
xlabel('Continuation steps N_c'); ylabel('||grad||');
title('Posterior mean: endpoint gradient vs steps'); legend([leg, {'low-fidelity optimizer'}]); grid on; drawnow;

figure;
loglog(t_run, err, '-o');
xlabel('Time (s)'); ylabel('Relative error');
title('Posterior mean: accuracy vs cost'); legend(leg); grid on; drawnow;

figure;
loglog(t_run, grad, '-o');
yline(res.lowfi_grad_norm, 'k--');
xlabel('Time (s)'); ylabel('||grad||');
title('Posterior mean: endpoint gradient vs cost'); legend([leg, {'low-fidelity optimizer'}]); grid on; drawnow;

disp('Mean Newton-polish plots done.');
