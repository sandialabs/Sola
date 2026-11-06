%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plots for Experiment 2c (Driver_CI_Curvature_Sweep), read from its saved
% results.  Bare-bones on purpose: plot/semilogy, labels, title, legend.

results_file = fullfile(fileparts(mfilename('fullpath')), 'ci-curvature-sweep-results.mat');
loaded = load(results_file, 'results');
res = loaded.results;

polish_Ms = res.settings.polish_Ms;
grad_target = res.settings.grad_target;

for ip = 1:length(polish_Ms)
    M = polish_Ms(ip);
    fx = res.fixed{ip};
    ad = res.adapt{ip};
    fam = {ad.family};
    both = strcmp(fam, 'both');
    curv_only = strcmp(fam, 'curv only');

    %% Median and tail gradient vs mean work
    fprintf('Plots (M = %d): endpoint gradient vs work\n', M);
    stats = {'median', 'p90'};
    titles = {'Median', '90th percentile'};
    for ss = 1:2
        figure;
        hold on;
        semilogy([fx.work_mean], [fx.(stats{ss})], 'ko-');
        semilogy([ad(both).work_mean], [ad(both).(stats{ss})], 's-');
        semilogy([ad(curv_only).work_mean], [ad(curv_only).(stats{ss})], '^-');
        Label_Points(ad(both), stats{ss});
        Label_Points(ad(curv_only), stats{ss});
        yline(grad_target, '--');
        set(gca, 'YScale', 'log');
        xlabel('Mean work per sample (steps + backtracks + Newton steps)');
        ylabel('Endpoint gradient norm');
        title(sprintf('%s gradient vs work (M = %d); labels = \\epsilon_{curv}', titles{ss}, M));
        legend('fixed N_c', sprintf('both (\\epsilon_{nl} = %g)', res.settings.eps_nl_both), 'curvature only', 'low-fidelity ||grad||', 'Location', 'best');
        drawnow;
    end

    %% Max gradient vs mean work
    figure;
    hold on;
    semilogy([fx.work_mean], [fx.max], 'ko-');
    semilogy([ad(both).work_mean], [ad(both).max], 's-');
    semilogy([ad(curv_only).work_mean], [ad(curv_only).max], '^-');
    Label_Points(ad(both), 'max');
    Label_Points(ad(curv_only), 'max');
    yline(grad_target, '--');
    set(gca, 'YScale', 'log');
    xlabel('Mean work per sample (steps + backtracks + Newton steps)');
    ylabel('Endpoint gradient norm');
    title(sprintf('Worst-sample gradient vs work (M = %d)', M));
    legend('fixed N_c', 'both', 'curvature only', 'low-fidelity ||grad||', 'Location', 'best');
    drawnow;
end

%% Does eta rank sample difficulty?  (1-step gradient vs eta)
disp('Plot: eta vs one-step gradient');
figure;
loglog(res.eta, res.one_step_grad, 'o');
xlabel('\eta (a priori curvature)');
ylabel('Gradient after one fixed step (unpolished)');
title(sprintf('Difficulty vs \\eta (Spearman %.2f)', res.pred_corr(1, 2)));
drawnow;

%% Predictor-defect coefficient vs eta and vs eta * ||beta_dot||
disp('Plot: predictor defect vs eta and eta*||beta_dot||');
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

function Label_Points(cfgs, stat)
    for kk = 1:length(cfgs)
        text(cfgs(kk).work_mean, cfgs(kk).(stat), sprintf('  %g', cfgs(kk).eps_curv), 'FontSize', 12);
    end
end
