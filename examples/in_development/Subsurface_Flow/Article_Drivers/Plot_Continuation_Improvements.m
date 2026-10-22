%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plots for Driver_Continuation_Improvements, read from its saved results
% file so figures can be regenerated without rerunning the solves:
%   Plot_Continuation_Improvements('continuation-improvements-results.mat', false)
%
% Text uses the default TeX interpreter.  Figures use the painters renderer and
% subplot/sgtitle instead of the newer tile layout, to stay clear of GPU/OpenGL paths.  Each figure is built in its own
% try/catch so a graphics failure in one does not prevent the others.

function Plot_Continuation_Improvements(results_file, export_figures)
    if nargin < 1 || isempty(results_file)
        results_file = fullfile(fileparts(mfilename('fullpath')), 'continuation-improvements-results.mat');
    end
    if nargin < 2
        export_figures = false;
    end
    out_dir = fileparts(results_file);
    loaded = load(results_file, 'results');
    res = loaded.results;

    figure_makers = {@Plot_Mean_To_Sample, 'improvements_mean_to_sample.pdf';};
    % figure_makers = {@Plot_Mean_To_Sample, 'improvements_mean_to_sample.pdf'; ...
    %                  @Plot_Sample_Bands, 'improvements_sample_bands.pdf'; ...
    %                  @Plot_Adaptive_Steps, 'improvements_adaptive_steps.pdf'};
    for ii = 1:size(figure_makers, 1)
        try
            fig = figure_makers{ii, 1}(res);
            drawnow;
            if export_figures
                exportgraphics(fig, fullfile(out_dir, figure_makers{ii, 2}));
            end
        catch ME
            warning('Plot_Continuation_Improvements:PlotFailed', 'Figure %d failed: %s', ii, ME.message);
        end
    end
end

function c = Colors()
    c.nominal = [0.85 0.33 0.10];
    c.mts = [0.00 0.45 0.74];
    c.adapt = [0.13 0.55 0.13];
    c.ref = [0.30 0.30 0.30];
end

function fig = Plot_Mean_To_Sample(res)
    c = Colors();
    steps = res.settings.step_counts;
    exp1 = res.exp1;
    num_samples = res.settings.num_samples;
    med = @(x) median(x, 'omitnan');
    hi = @(x) prctile(x(isfinite(x)), 90);
    nom_med = arrayfun(@(e) med(e.nominal.ref_err), exp1);
    nom_hi = arrayfun(@(e) hi(e.nominal.ref_err), exp1);
    mts_med = arrayfun(@(e) med(e.mts.ref_err), exp1);
    mts_hi = arrayfun(@(e) hi(e.mts.ref_err), exp1);
    nom_grad = arrayfun(@(e) med(e.nominal.grad_norm), exp1);
    mts_grad = arrayfun(@(e) med(e.mts.grad_norm), exp1);
    nom_time = arrayfun(@(e) e.nominal.time, exp1) / num_samples;
    % Mean-to-sample cost includes its one-time mean solve, amortized over samples.
    mts_time = arrayfun(@(e) e.mts.time, exp1) / num_samples + res.mean_start.time / num_samples;

    fig = figure('Name', 'Mean-to-sample vs nominal-to-sample', 'Position', [100 100 1300 900], 'Renderer', 'painters');
    nrow = 2; ncol = 2; tile = 0;
    sgtitle('Mean-to-sample continuation reaches a given accuracy in far fewer steps');

    tile = tile + 1; subplot(nrow, ncol, tile);
    semilogy(steps, nom_med, 'o-', 'Color', c.nominal, 'DisplayName', 'nominal-to-sample (median)'); hold on;
    semilogy(steps, nom_hi, 'o--', 'Color', c.nominal, 'DisplayName', 'nominal-to-sample (90th pct)');
    semilogy(steps, mts_med, 's-', 'Color', c.mts, 'DisplayName', 'mean-to-sample (median)');
    semilogy(steps, mts_hi, 's--', 'Color', c.mts, 'DisplayName', 'mean-to-sample (90th pct)');
    grid on;
    xlabel('sample continuation steps N_c');
    ylabel('relative error vs. reference');
    title('Accuracy at equal steps');
    legend('Location', 'southwest', 'FontSize', 12);

    tile = tile + 1; subplot(nrow, ncol, tile);
    semilogy(steps, nom_grad, 'o-', 'Color', c.nominal, 'DisplayName', 'nominal-to-sample'); hold on;
    semilogy(steps, mts_grad, 's-', 'Color', c.mts, 'DisplayName', 'mean-to-sample');
    grid on;
    xlabel('sample continuation steps N_c');
    ylabel('median endpoint ||\nabla_\beta J||');
    title('Optimality of the sample (no reference needed)');
    legend('Location', 'southwest', 'FontSize', 12);

    tile = tile + 1; subplot(nrow, ncol, tile);
    semilogy(nom_time, nom_med, 'o-', 'Color', c.nominal, 'DisplayName', 'nominal-to-sample'); hold on;
    semilogy(mts_time, mts_med, 's-', 'Color', c.mts, 'DisplayName', 'mean-to-sample (incl. mean solve)');
    grid on;
    xlabel('wall time per sample (s)');
    ylabel('median relative error');
    title('Accuracy at equal cost');
    legend('Location', 'southwest', 'FontSize', 12);

    tile = tile + 1; subplot(nrow, ncol, tile);
    [~, jshow] = min(abs(steps - res.settings.show_N_c));
    e_nom = exp1(jshow).nominal.ref_err;
    e_mts = exp1(jshow).mts.ref_err;
    Plot_Paired_Scatter(e_nom, e_mts, c.mts, 'nominal-to-sample error', 'mean-to-sample error');
    title(sprintf('Per sample, N_c = %d (below line = mean-to-sample wins)', steps(jshow)));
end

function fig = Plot_Sample_Bands(res)
    c = Colors();
    ptile = 90;
    steps = res.settings.step_counts;
    [~, jshow] = min(abs(steps - res.settings.show_N_c));
    x_plot = 1:numel(res.z_opt);
    ref_band = Band(res.sample_ref.z, ptile);
    panels = {res.sample_ref.z, sprintf('Reference (mean-to-sample, N_c = %d)', res.settings.reference_steps), c.ref; ...
              res.exp1(jshow).nominal.z, 'Nominal-to-sample', c.nominal; ...
              res.exp1(jshow).mts.z, 'Mean-to-sample', c.mts};

    fig = figure('Name', 'Sample bands at small step budget', 'Position', [100 100 1500 480], 'Renderer', 'painters');
    nrow = 1; ncol = 3; tile = 0;
    sgtitle(sprintf('Posterior samples with N_c = %d steps vs. reference', steps(jshow)));
    for pp = 1:3
        tile = tile + 1; subplot(nrow, ncol, tile); hold on;
        lh = Band(panels{pp, 1}, ptile);
        fill([x_plot, fliplr(x_plot)], [lh(:, 1).', fliplr(lh(:, 2).')], panels{pp, 3}, ...
             'FaceAlpha', 0.4, 'EdgeColor', 'none', 'DisplayName', 'samples');
        plot(x_plot, ref_band(:, 1), ':', 'Color', c.ref, 'LineWidth', 2, 'DisplayName', 'reference band');
        plot(x_plot, ref_band(:, 2), ':', 'Color', c.ref, 'LineWidth', 2, 'HandleVisibility', 'off');
        plot(x_plot, res.z_hifi, 'r-', 'LineWidth', 2, 'DisplayName', 'HiFi');
        xlabel('Control index');
        title(panels{pp, 2});
        if pp == 1, ylabel('Injection/production rate'); end
        legend('Location', 'best', 'FontSize', 11);
    end
end

function fig = Plot_Adaptive_Steps(res)
    c = Colors();
    num_samples = res.settings.num_samples;
    mean_fixed = res.mean_fixed;
    mean_adapt = res.mean_adapt;
    sample_fixed = res.sample_fixed;
    adapt = res.sample_adapt{res.settings.show_idx};
    hi = @(x) prctile(x(isfinite(x)), 90);

    fig = figure('Name', 'Adaptive time steps', 'Position', [100 100 1300 900], 'Renderer', 'painters');
    nrow = 2; ncol = 2; tile = 0;
    sgtitle('A priori time-step heuristic: more accuracy per unit of work');

    tile = tile + 1; subplot(nrow, ncol, tile);
    semilogy([mean_fixed.work], [mean_fixed.ref_err], 'o-', 'Color', c.mts, 'DisplayName', 'fixed uniform steps'); hold on;
    semilogy([mean_adapt.work], [mean_adapt.ref_err], 's-', 'Color', c.adapt, 'MarkerFaceColor', c.adapt, ...
             'MarkerSize', 12, 'DisplayName', 'adaptive (varying \epsilon_{curv})');
    grid on;
    xlabel('work (accepted steps + backtracks)');
    ylabel('relative error vs. reference');
    title('Posterior mean');
    legend('Location', 'northeast', 'FontSize', 12);

    tile = tile + 1; subplot(nrow, ncol, tile);
    fixed_work = res.settings.sample_fixed_sweep;
    fixed_med = cellfun(@(s) median(s.ref_err, 'omitnan'), sample_fixed);
    fixed_hi = cellfun(@(s) hi(s.ref_err), sample_fixed);
    semilogy(fixed_work, fixed_med, 'o-', 'Color', c.mts, 'DisplayName', 'fixed (median)'); hold on;
    semilogy(fixed_work, fixed_hi, 'o--', 'Color', c.mts, 'DisplayName', 'fixed (90th pct)');
    all_adapt = res.sample_adapt;
    adapt_work = cellfun(@(a) mean(a.work, 'omitnan'), all_adapt);
    semilogy(adapt_work, cellfun(@(a) median(a.ref_err, 'omitnan'), all_adapt), 's-', 'Color', c.adapt, ...
             'MarkerFaceColor', c.adapt, 'MarkerSize', 12, 'DisplayName', 'adaptive (median)');
    semilogy(adapt_work, cellfun(@(a) hi(a.ref_err), all_adapt), 's--', 'Color', c.adapt, ...
             'MarkerSize', 12, 'DisplayName', 'adaptive (90th pct)');
    grid on;
    xlabel('mean work per sample');
    ylabel('relative error vs. reference');
    title('Mean-to-sample, across samples');
    legend('Location', 'northeast', 'FontSize', 12);

    tile = tile + 1; subplot(nrow, ncol, tile);
    bar(1:num_samples, [adapt.steps(:), adapt.backtracks(:)], 'stacked');
    grid on;
    xlabel('sample index');
    ylabel('steps used');
    title(sprintf('Per-sample work (\\epsilon_{curv} = %g)', res.settings.sample_epsilon_curv(res.settings.show_idx)));
    legend({'accepted steps', 'backtracks'}, 'Location', 'best', 'FontSize', 12);

    tile = tile + 1; subplot(nrow, ncol, tile);
    fixed_p90 = cellfun(@(s) hi(s.grad_norm), sample_fixed);
    fixed_max = cellfun(@(s) max(s.grad_norm, [], 'omitnan'), sample_fixed);
    adapt_p90 = cellfun(@(a) hi(a.grad_norm), all_adapt);
    adapt_max = cellfun(@(a) max(a.grad_norm, [], 'omitnan'), all_adapt);
    semilogy(fixed_work, fixed_p90, 'o-', 'Color', c.mts, 'DisplayName', 'fixed (90th pct)'); hold on;
    semilogy(fixed_work, fixed_max, 'o--', 'Color', c.mts, 'DisplayName', 'fixed (max)');
    semilogy(adapt_work, adapt_p90, 's-', 'Color', c.adapt, 'DisplayName', 'adaptive (90th pct)');
    semilogy(adapt_work, adapt_max, 's--', 'Color', c.adapt, 'DisplayName', 'adaptive (max)');
    grid on;
    xlabel('mean work per sample');
    ylabel('endpoint gradient norm');
    title('Worst-case samples');
    legend('Location', 'northeast', 'FontSize', 12);
end

function Plot_Paired_Scatter(x, y, color, x_label, y_label)
    % Log-log scatter against the y = x line, annotated with the win count.
    ok = isfinite(x) & isfinite(y) & x > 0 & y > 0;
    hold on;
    if any(ok)
        lims = [min([x(ok); y(ok)]) / 2, max([x(ok); y(ok)]) * 2];
        loglog(lims, lims, 'k--', 'LineWidth', 1.5);
        loglog(x(ok), y(ok), 'o', 'Color', color, 'MarkerFaceColor', color, 'MarkerSize', 9);
        xlim(lims); ylim(lims);
    end
    set(gca, 'XScale', 'log', 'YScale', 'log'); grid on;
    xlabel(x_label);
    ylabel(y_label);
    text(0.05, 0.92, sprintf('wins on %d / %d samples', sum(y(ok) < x(ok)), sum(ok)), ...
         'Units', 'normalized', 'FontSize', 14);
end

function lh = Band(samples, ptile)
    lh = [prctile(samples, 50 - ptile / 2, 2), prctile(samples, 50 + ptile / 2, 2)];
end
