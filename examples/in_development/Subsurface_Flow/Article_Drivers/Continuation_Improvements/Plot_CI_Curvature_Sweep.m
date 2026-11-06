%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plots for Experiment 2c (Driver_CI_Curvature_Sweep), read from its saved
% results.  With Newton polishing that stops on a gradient tolerance, the
% endpoint gradient is small by construction, so the plots compare the cost to
% reach the tolerance and the fraction of samples that reach it.

results_file = fullfile(fileparts(mfilename('fullpath')), 'ci-curvature-sweep-results.mat');
loaded = load(results_file, 'results');
res = loaded.results;

cost_label = sprintf('Cost per sample (continuation steps + %.2g \\times Newton steps)', res.settings.newton_cost);
header = sprintf('Mean-to-sample, Newton until ||g|| \\leq %.2g', res.settings.newton_gradient_tol);

%% 1-2. Cost to reach the tolerance: fixed N_c vs adaptive settings
disp('Plots 1-2: cost to reach the Newton tolerance');
fixed_summaries = cellfun(@(f) f.summary, res.fixed, 'UniformOutput', false);
adapt_summaries = cellfun(@(a) a.summary, res.adapt, 'UniformOutput', false);
summaries = [fixed_summaries{:}, adapt_summaries{:}];
Plot_CI_Cost_Panels(summaries, header, cost_label);

%% 3. Does eta rank sample difficulty?  (Newton steps needed after one step vs eta)
disp('Plot 3: eta vs difficulty');
figure;
semilogx(res.eta, res.difficulty, 'o');
xlabel('\eta (a priori curvature)');
ylabel('Newton steps needed after one continuation step');
title(sprintf('Difficulty vs \\eta (Spearman %.2f)', res.pred_corr(1, 2)));
drawnow;

%% 4-5. Predictor-defect coefficient vs eta and vs eta * ||beta_dot||
disp('Plots 4-5: predictor defect vs eta and eta*||beta_dot||');
c_pred = res.pred.c(:, end);
figure;
loglog(res.eta, c_pred, 'o');
xlabel('\eta');
ylabel('Predictor defect coefficient c = ||g(pred) - g(start)|| / dt^2');
title(sprintf('Predictor defect vs \\eta (Spearman %.2f)', res.pred_corr(1, 1)));
drawnow;

figure;
loglog(res.eta .* res.pred.speed, c_pred, 'o');
xlabel('\eta \times ||\beta_{dot}||');
ylabel('Predictor defect coefficient c = ||g(pred) - g(start)|| / dt^2');
title(sprintf('Predictor defect vs \\eta||\\beta_{dot}|| (Spearman %.2f)', res.pred_corr(2, 1)));
drawnow;

disp('Curvature-sweep plots done.');
