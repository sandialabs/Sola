%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot for Experiment 2a (Driver_CI_Adaptive_Mean), read from its saved results.

results_file = fullfile(fileparts(mfilename('fullpath')), 'ci-adaptive-mean-results.mat');
loaded = load(results_file, 'results');
res = loaded.results;

%% 4. Adaptive steps, posterior mean
disp('Plot 4: adaptive steps, mean');
figure;
semilogy([res.mean_fixed.work], [res.mean_fixed.grad_norm], 'o-', ...
         [res.mean_adapt.work], [res.mean_adapt.grad_norm], 's-');
xlabel('Work (steps + backtracks + Newton steps)');
ylabel('Endpoint gradient norm');
title('Posterior mean: fixed vs adaptive');
legend('fixed', 'adaptive');
drawnow;

disp('Adaptive-mean plots done.');
