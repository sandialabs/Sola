%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plots for Experiment 2a (Driver_CI_Adaptive_Mean), read from its saved results.
% Bare-bones on purpose: plot/semilogy, labels, title, legend.

results_file = fullfile(fileparts(mfilename('fullpath')), 'ci-adaptive-mean-results.mat');
loaded = load(results_file, 'results');
res = loaded.results;

runs = [res.mean_fixed, res.mean_backtrack, res.mean_adapt];
labels = {runs.label};
n = length(runs);

%% 1. Cost to reach the Newton tolerance: continuation vs Newton
disp('Plot 1: cost of continuation and Newton steps');
cont = [runs.steps] + [runs.backtracks];
newton_cost_part = [runs.cost] - cont;
figure;
bar(1:n, [cont(:), newton_cost_part(:)], 'stacked');
set(gca, 'XTick', 1:n, 'XTickLabel', labels, 'TickLabelInterpreter', 'none');
xtickangle(45);
ylabel(sprintf('Cost (continuation steps + %.2g \\times Newton steps)', res.settings.newton_cost));
title(sprintf('Posterior mean, Newton until ||g|| \\leq %.2g', res.settings.newton_gradient_tol));
legend('continuation steps', 'Newton steps (weighted)', 'Location', 'best');
drawnow;

%% 2. Newton convergence trajectories
disp('Plot 2: gradient norm vs Newton step');
figure;
hold on;
for k = 1:n
    g = runs(k).newton_grad_norms;
    semilogy(0:(length(g) - 1), g, 'o-');
end
yline(res.settings.newton_gradient_tol, '--');
set(gca, 'YScale', 'log');
xlabel('Newton steps after the continuation');
ylabel('Gradient norm');
title('Posterior mean: Newton convergence from each continuation');
legend([labels, {'tolerance'}], 'Location', 'best', 'Interpreter', 'none');
drawnow;

disp('Adaptive-mean plots done.');
