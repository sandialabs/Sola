%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plots for Experiment 2d (Driver_CI_Polish_Sweep), read from its saved
% results.  Bare-bones on purpose: plot/semilogy, labels, title, legend.

results_file = fullfile(fileparts(mfilename('fullpath')), 'ci-polish-sweep-results.mat');
loaded = load(results_file, 'results');
res = loaded.results;

cfgs = res.cfgs;
num_cfg = length(cfgs);
grad_target = res.settings.grad_target;
names = {cfgs.name};

%% 1. Gradient vs number of Newton steps
disp('Plot 1: gradient vs Newton steps');
stats = {'median', 'p90'};
titles = {'Median', '90th percentile'};
for ss = 1:2
    figure;
    hold on;
    for jj = 1:num_cfg
        semilogy([cfgs(jj).fixedM.m], [cfgs(jj).fixedM.(stats{ss})], 'o-');
    end
    yline(grad_target, '--');
    set(gca, 'YScale', 'log');
    xlabel('Newton polishing steps M');
    ylabel('Endpoint gradient norm');
    title([titles{ss} ' gradient vs Newton steps']);
    legend([names, {'low-fidelity ||grad||'}], 'Location', 'best');
    drawnow;
end

%% 2. Gradient vs cost; fixed M (circles, labeled by M) and stopping rule (stars)
disp('Plot 2: gradient vs cost');
for ss = 1:2
    figure;
    hold on;
    hs = gobjects(1, num_cfg);
    for jj = 1:num_cfg
        h = semilogy([cfgs(jj).fixedM.cost], [cfgs(jj).fixedM.(stats{ss})], 'o-');
        hs(jj) = h;
        semilogy([cfgs(jj).stopM.cost], [cfgs(jj).stopM.(stats{ss})], 'p', 'Color', h.Color, 'MarkerFaceColor', h.Color);
        for k = 1:length(cfgs(jj).fixedM)
            text(cfgs(jj).fixedM(k).cost, cfgs(jj).fixedM(k).(stats{ss}), sprintf('  %d', cfgs(jj).fixedM(k).m), 'FontSize', 12);
        end
    end
    yline(grad_target, '--');
    set(gca, 'YScale', 'log');
    xlabel(sprintf('Mean cost (steps + %.2g Newton steps)', res.settings.newton_cost));
    ylabel('Endpoint gradient norm');
    title([titles{ss} ' gradient vs cost (labels = M, stars = stopping rule)']);
    legend(hs, names, 'Location', 'best');
    drawnow;
end

disp('Polish-sweep plots done.');
