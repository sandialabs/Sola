%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Plot_CI_Cost_Panels(summaries, header, cost_label)
    % Two figures comparing continuation settings by the work needed to reach
    % the Newton gradient tolerance (see Summarize_CI_Cost):
    %   1. stacked bars of the mean continuation and Newton cost per sample,
    %      with the 90th-percentile and max sample cost, and the converged
    %      fraction written above each bar;
    %   2. empirical CDF of the per-sample cost (tails = hard samples).
    n = length(summaries);
    labels = {summaries.label};

    mean_cont = [summaries.mean_cont];
    mean_newton_cost = [summaries.mean_cost] - mean_cont;
    p90 = [summaries.p90_cost];
    mx = [summaries.max_cost];
    conv = [summaries.converged_frac];

    figure;
    bar(1:n, [mean_cont(:), mean_newton_cost(:)], 'stacked');
    hold on;
    plot(1:n, p90, 'kv', 'MarkerFaceColor', 'k');
    plot(1:n, mx, 'k^', 'MarkerFaceColor', 'k');
    for k = 1:n
        text(k, mx(k), sprintf('%.0f%%', 100 * conv(k)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 12);
    end
    set(gca, 'XTick', 1:n, 'XTickLabel', labels, 'TickLabelInterpreter', 'none');
    xtickangle(45);
    ylabel(cost_label);
    title([header ' (labels: % of samples reaching the tolerance)']);
    legend('continuation steps', 'Newton steps (weighted)', '90th percentile', 'max', 'Location', 'best');
    drawnow;

    figure;
    hold on;
    legend_labels = cell(1, n);
    for k = 1:n
        c = sort(summaries(k).cost(isfinite(summaries(k).cost)));
        stairs([c; c(end)], [(1:length(c))' / length(c); 1], '-');
        legend_labels{k} = sprintf('%s (%.0f%% conv.)', labels{k}, 100 * conv(k));
    end
    xlabel(cost_label);
    ylabel('Fraction of samples with cost \leq x');
    title([header ': per-sample cost distribution']);
    legend(legend_labels, 'Location', 'southeast', 'Interpreter', 'none');
    drawnow;
end
