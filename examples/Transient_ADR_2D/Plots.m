clear;
close all;
clc;
run('../../src/Set_Paths');

%% Plot configuration.
set(0, 'DefaultAxesFontName', 'Times New Roman');
set(0, 'DefaultAxesFontWeight', 'normal');
set(0, 'DefaultAxesFontSize', 20);
set(0, 'DefaultTextFontName', 'Times New Roman');
set(0, 'DefaultTextFontWeight', 'normal');
set(0, 'DefaultTextFontSize', 32);
cmap = viridis;

%% Load solver and other data.
% model = Transient_ADR_2D.model_fromfile(meshfile);
load('OpInf_Training_Data.mat', 'Z_train');
load('OptimizationSolution.mat', 't', 'solver', 'basis1', 'basis2', 'obj_hifi', 'Q_rom', 'Y_hifi', 'rs', 'opt');
load('fem_matrices.mat', 'mass_matrix');
load('MD_Results.mat');

%% Finite element mesh.
fig = plotfield(solver, "mesh", 0, true, "Spatial mesh and injection locations");
print(fig, 'figures/adr_mesh.pdf', '-dpdf', '-r300', '-loose');
close(fig);

%% Protection zone.
fig = plotfield(solver, obj_hifi.target_weight, viridis, true, "Protection zone weights \psi");
print(fig, 'figures/adr_target.png', '-dpng', '-r300', '-loose');
close(fig);

%% Initial condition.
u10 = 50 * exp(-50 .* sum(([solver.x, solver.y]' - solver.init_center).^2, 1))';
fig = plotfield(solver, u10, viridis, true, "Initial condition for contaminant");
print(fig, 'figures/adr_initial.png', '-dpng', '-r300', '-loose');
close(fig);

%% POD singular value decay for each state variable.
svdvals = [basis1.singular_values; basis2.singular_values];
j = 0:size(svdvals, 2);
resenergy = 1 - (cumsum(svdvals.^2, 2) ./ sum(svdvals.^2, 2));
resenergy = [1, resenergy(1, :); 1, resenergy(2, :)];

fig = figure;
ax = subplot(1, 1, 1);
for i = [1, 2]
    semilogy(ax, j, resenergy(i, :), '.-', LineWidth = 1.5, MarkerSize = 15);
    hold on;
end
xline(ax, rs(1), 'Color', ax.ColorOrder(1, :));
xline(ax, rs(2), 'Color', ax.ColorOrder(2, :));
yline(ax, 1e-5);
xlim(ax, [0, 50]);
ylim(ax, [1e-8, 1e0]);
xlabel(ax, 'Singular value index');
ylabel(ax, 'Residual energy');
text(ax, rs(1) - 0.5, 2e-8, '$r_1 = 24$', 'Interpreter', 'latex', ...
     'HorizontalAlignment', 'right', 'Color', ax.ColorOrder(1, :), ...
     'FontSize', 16);
text(ax, rs(2) + 0.5, 2e-8, '$r_2 = 30$', 'Interpreter', 'latex', ...
     'HorizontalAlignment', 'left', 'Color', ax.ColorOrder(2, :), ...
     'FontSize', 16);
set(fig, 'Position', [175, 300, 560, 330]);
print(fig, 'figures/adr_svdvals.pdf', '-dpdf', '-r300', '-loose');
close(fig);

%% Example training controls
Q1 = reshape(Z_train(:, 1), solver.n_q, []);
Q1 = [zeros(solver.n_q, 1), Q1];
nlines = size(Q1, 1);
fig = figure;
ax = subplot(1, 1, 1);
% colors = lines(nlines);
colors = distinguishable_colors(nlines);
hold on;
for i = 1:nlines
    plot(ax, t, Q1(i, :), 'Color', colors(i, :), 'LineWidth', 1.5);
end
xlabel(ax, '$t$', 'Interpreter', 'latex');
ylabel(ax, '$q_i^{(\ell)}(t)$', 'Interpreter', 'latex');
ylim(ax, [0, 15]);
set(fig, 'Position', [175, 300, 560, 330]);
set(title('Training Controls'), 'FontWeight', 'normal');
print(fig, 'figures/adr_traincontrols.pdf', '-dpdf', '-r300', '-loose');
close(fig);

%% Dominant basis vectors for each variable.
% [fig, ax] = plotfield(solver, basis1.singular_vectors(:, 1), viridis, true, "First basis function for contaminant");
% print(fig, 'figures/adr_basis1.png', '-dpng', '-r300', '-loose');
% close(fig);
%
% [fig, ax] = plotfield(solver, basis2.singular_vectors(:, 1), "parula", true, "First basis function for neutralizer");
% print(fig, 'figures/adr_basis2.png', '-dpng', '-r300', '-loose');
% close(fig);

%% Grid of basis vectors over the spatial domain.
nvecs = 4;
svdvecs1 = basis1.singular_vectors(:, 1:nvecs);
svdvecs2 = basis2.singular_vectors(:, 1:nvecs);
lim1 = max(abs([min(svdvecs1(:)), max(svdvecs2(:))]));
lim2 = max(abs([min(svdvecs2(:)), max(svdvecs2(:))]));
for j = 1:nvecs
    ytextA = '';
    ytextB = '';
    if j == 1
        ytextA = 'Contaminant';
        ytextB = 'Neutralizer';
    end
    [figA, axA] = plotfield(solver, svdvecs1(:, j), "viridis", j == 4, ['$j=' num2str(j) '$'], ytextA, [-lim1, lim1]);
    print(figA, ['figures/adr_basis1-', num2str(j), '.png'], '-dpng', '-r300', '-loose');
    close(figA);

    [figB, axB] = plotfield(solver, svdvecs2(:, j), "parula", j == 4, '', ytextB, [-lim2, lim2]);
    print(figB, ['figures/adr_basis2-', num2str(j), '.png'], '-dpng', '-r300', '-loose');
    close(figB);
end

%% ROMCO control solution as a function of time.
fig = figure;
ax = subplot(1, 1, 1);
hold on;
for i = 1:nlines
    plot(ax, t(2:end), abs(Q_rom(i, :)), 'Color', colors(i, :), 'LineWidth', 1.5);
end
xlabel(ax, '$t$', 'Interpreter', 'latex');
ylabel(ax, '$q_i^{(\ell)}(t)$', 'Interpreter', 'latex');
ylim(ax, [0, 15]);
set(fig, 'Position', [175, 300, 560, 330]);
set(title('Optimized Controls'), 'FontWeight', 'normal');
print(fig, 'figures/adr_optcontrols.pdf', '-dpdf', '-r300', '-loose');
close(fig);

%% FOM solution with ROMCO controls as a function of time (including init).
nsnaps = 4;
Y1 = Y_hifi(:, 1, :);
Y2 = Y_hifi(:, 2, :);
lim1 = 20; % max(Y1(:)) / 2;
lim2 = 20; % max(Y2(:)) / 1.5;
limmin = 1e-1;
indices = round(linspace(1, length(t), nsnaps + 1));
indices = indices(2:end);
for j = 1:nsnaps
    ytextA = '';
    ytextB = '';
    if j == 1
        ytextA = 'Contaminant';
        ytextB = 'Neutralizer';
    end
    [figA, axA] = plotfield(solver, Y1(:, indices(j)), "viridis", j == nsnaps, ['$t=' num2str(t(indices(j))) '$'], ytextA, [limmin, lim1], false, true);
    print(figA, ['figures/adr_romcofom1-', num2str(j), '.png'], '-dpng', '-r300', '-loose');
    close(figA);

    [figB, axB] = plotfield(solver, Y2(:, indices(j)), "parula", j == nsnaps, '', ytextB, [limmin, lim2], false, true);
    print(figB, ['figures/adr_romcofom2-', num2str(j), '.png'], '-dpng', '-r300', '-loose');
    close(figB);
end

%%
solver.vel_params = load('OptimizationSolution.mat', 'vel_params_rom').vel_params_rom;
solver.Plot_Velocity_Field(20);
xticks([]);
yticks([]);
set(title('Training Wind Field'), 'FontWeight', 'normal');
print(figure(1), 'figures/training_wind.pdf', '-dpdf', '-r300', '-loose');
close;

solver.vel_params = load('OptimizationSolution.mat', 'vel_params_hifi').vel_params_hifi;
solver.Plot_Velocity_Field(20);
xticks([]);
yticks([]);
set(title('Testing Wind Field'), 'FontWeight', 'normal');
print(figure(1), 'figures/testing_wind.pdf', '-dpdf', '-r300', '-loose');
close;

%% MD updated control solution as a function of time.
n_t = length(t);
n_q = length(z_update_mean) / (n_t - 1);
Q_update_mean = reshape(z_update_mean, n_q, n_t - 1);
[~, I] = sort(vecnorm((abs(Q_update_mean) - abs(Q_rom))'), 'descend');
fig = figure;
ax = subplot(1, 1, 1);
hold on;
for i = I(1:3)
    plot(ax, t(2:end), abs(Q_update_mean(i, :)), '--', 'Color', colors(i, :), 'LineWidth', 1.5);
    plot(ax, t(2:end), abs(Q_rom(i, :)), 'Color', colors(i, :), 'LineWidth', 1.5);
end
xlabel(ax, '$t$', 'Interpreter', 'latex');
ylabel(ax, '$q_i^{(\ell)}(t)$', 'Interpreter', 'latex');
set(fig, 'Position', [175, 300, 560, 330]);
print(fig, 'figures/adr_updatedoptcontrols.pdf', '-dpdf', '-r300', '-loose');
close(fig);

% i = 11;
% num_post_samples = size(z_update_samples,2);
% Q_update_samples = zeros(n_q,n_t-1,num_post_samples);
% for k = 1:num_post_samples
%     Q_update_samples(:,:,k) = reshape(z_update_samples(:,k), n_q, n_t - 1);
% end
% fig = figure;
% ax = subplot(1, 1, 1);
% hold on;
% plot(ax, t(2:end), abs(Q_update_mean(i, :)), '--', 'Color', colors(i, :), 'LineWidth', 1.5);
% plot(ax, t(2:end), abs(Q_rom(i, :)), 'Color', colors(i, :), 'LineWidth', 1.5);
% for k = 1:num_post_samples
%    plot(ax,t(2:end),abs(Q_update_samples(i,:,k)),'Color',[.9,.9,.9],'LineWidth',1.5);
% end
% plot(ax, t(2:end), abs(Q_update_mean(i, :)), '--', 'Color', colors(i, :), 'LineWidth', 1.5);
% plot(ax, t(2:end), abs(Q_rom(i, :)), 'Color', colors(i, :), 'LineWidth', 1.5);
% xlabel(ax, '$t$', 'Interpreter', 'latex');
% ylabel(ax, '$q_{11}^{(\ell)}(t)$', 'Interpreter', 'latex');

%% FOM solution with MD updated controls as a function of time (including init).
nsnaps = 5;
Y1 = Y_update_mean(:, 1, :);
Y2 = Y_update_mean(:, 2, :);
lim1 = 20; % max(Y1(:)) / 2;
lim2 = 20; % max(Y2(:)) / 1.5;
limmin = 1e-1;
indices = round(linspace(1, length(t), nsnaps));
for j = 1:nsnaps
    ytextA = '';
    ytextB = '';
    if j == 1
        ytextA = 'Contaminant';
        ytextB = 'Neutralizer';
    end
    [figA, axA] = plotfield(solver, Y1(:, indices(j)), "viridis", j == nsnaps, ['$t=' num2str(t(indices(j))) '$'], ytextA, [limmin, lim1], false);
    print(figA, ['figures/adr_updatedfom1-', num2str(j), '.png'], '-dpng', '-r300', '-loose');
    close(figA);

    [figB, axB] = plotfield(solver, Y2(:, indices(j)), "parula", j == nsnaps, '', ytextB, [limmin, lim2], false);
    print(figB, ['figures/adr_updatedfom2-', num2str(j), '.png'], '-dpng', '-r300', '-loose');
    close(figB);
end

%%
val_hifi = zeros(n_t, 1);
val_update = zeros(n_t, 1);
u_hifi = Y_hifi(:, 1, :);
u_update = Y_update_mean(:, 1, :);
for k = 1:n_t
    val_hifi(k) = opt.obj.fullobj.g(u_hifi(:, k), t(k));
    val_update(k) = opt.obj.fullobj.g(u_update(:, k), t(k));
end

fig = figure;
ax = subplot(1, 1, 1);
hold on;
plot(ax, t, val_hifi, 'Color', colors(1, :), 'LineWidth', 1.5);
plot(ax, t, val_update, '--', 'Color', colors(2, :), 'LineWidth', 1.5);
xlabel(ax, '$t$', 'Interpreter', 'latex');
ylabel(ax, '$\| \mathbf{u}_1(t)*\mathbf{p} \|_{\mathbf{M}}^2$', 'Interpreter', 'latex');
legend({'ROMCO Controls', 'Posterior Mean Controls'}, 'Location', 'northwest');
set(fig, 'Position', [175, 300, 560, 330]);
print(fig, 'figures/contaminant_mass.pdf', '-dpdf', '-r300', '-loose');
close(fig);

%%
fig = figure;
ax = subplot(1, 1, 1);
hold on;
histogram(ax, obj_update_samples, 10, 'Normalization', 'probability');
h1 = plot(ax, [obj_lofi, obj_lofi], [0, ax.YLim(2)], 'LineWidth', 3);
h2 = plot(ax, [obj_update_mean, obj_update_mean], [0, ax.YLim(2)], 'LineWidth', 3);
xlabel(ax, 'High-fidelity objective value', 'Interpreter', 'latex');
ylabel(ax, 'Probability', 'Interpreter', 'latex');
legend([h1, h2], {'ROMCO Controls', 'Posterior Mean Controls'});
print(fig, 'figures/hifi_obj_dist.pdf', '-dpdf', '-r300', '-loose');
close(fig);

%% Postprocess saved figures.
[os, ~] = computer;
if strcmp(os, 'MACA64') || strcmp(os, 'MACI64')
    % Postprocessing (on MacOS):
    %  - Download TeX from https://www.tug.org/mactex/mactex-download.html
    %  - Install brew from https://brew.sh/
    %  - Check `which pdfcrop` and `which magick` give the paths below.
    setenv('PATH', [getenv('PATH') ':/Library/TeX/texbin' ':/opt/homebrew/bin' ':/usr/local/opt/ghostscript/bin' ':/usr/local/Cellar/imagemagick/7.1.1-39/bin/']);
    !for fname in figures/*.pdf; do pdfcrop --margin "1" "${fname}" "${fname}"; done
    !for fname in figures/*.png; do magick "${fname}" -trim +repage -bordercolor White -border 10x20 "${fname}"; done
    combinepngs('basis');
    combinepngs('romcofom');
    disp('ALL FIGURES PROCESSED AND SAVED');
else
    disp('Postprocessing only defined for MacOS');
end

%% Helper functions

function [fig, ax] = plotfield(solver, data, cmp, cb, titletext, ylabeltext, limits, logscale, setticks)
    % data - XY data to plot
    % cmp - colormap (viridis or 'parula')
    % cb - whether to include a colorbar (true or false)
    arguments
        solver
        data                % XY data to plot
        cmp = viridis       % colormap (viridis or 'parula')
        cb = true           % whether to include a colorbar (true or false)
        titletext = ''      % text for the title
        ylabeltext = ''     % text for the ylabel
        limits = []         % limits for the data
        logscale = false    % if True, logarithmically scale the colors.
        setticks = false    % if True, set the ticks of the colorbar manually.
    end

    fig = figure;
    ax = subplot(1, 1, 1);

    % Plot the data, add colorbar, etc.
    if isstring(data) && data == "mesh"
        pdemesh(solver.model);
        hold on;
        scatter(ax, solver.control_nodes(1, :), solver.control_nodes(2, :), 36, "black", "filled", "o");
    else
        fill(ax, [0 1.2 1.2 0], [0 0 1.2 1.2], [0.6 0.6 0.6], 'EdgeColor', 'none');
        ax.Layer = 'top';
        hold on;
        pdeplot(solver.model.Mesh, XYData = data, ColorMap = cmp, ColorBar = "off");
        if cb
            % Colorbar.
            cbar = colorbar(ax);
            colorbarPos = cbar.Position;
            colorbarPos(1) = colorbarPos(1) - 0.01;
            cbar.Position = colorbarPos;
            if setticks
                cbar.Ticks = 2:2:18;
            end
        end
        if ~isempty(limits)
            % Function call depends on Matlab version
            % clim(ax, limits);
            caxis(ax, limits);
        end
        if logscale
            set(ax, 'ColorScale', 'log');
        end
    end

    % Axis limits, etc.
    format_ax(solver, ax);

    % Title / labels
    if ~isempty(titletext)
        if contains(titletext, "$")
            title(ax, titletext, 'Interpreter', 'latex');
        else
            set(title(ax, titletext), 'FontWeight', 'normal');
        end
    end

    if ~isempty(ylabeltext)
        ylabelHandle = ylabel(ax, ylabeltext);
        ypos = get(ylabelHandle, 'Position');
        text(ypos(1) + 0.15, ypos(2), ylabeltext, 'Rotation', 90, ...
             'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
             'FontSize', 32);
    end
end

function format_ax(solver, ax)
    x = solver.x(:);
    y = solver.y(:);
    xlim(ax, [min(x) max(y)]);
    ylim(ax, [min(y) max(y)]);
    axis(ax, 'equal');
    box(ax, 'off');
    set(ax, 'XTick', [], 'YTick', [], 'XTickLabel', [], 'YTickLabel', []);
    set(ax, 'Color', 'none', 'XColor', 'none', 'YColor', 'none');
    set(ax, 'Position', [0.1300 0.1100 0.7750 0.8150]);

end

function combinepngs(label)
    c1 = ['magick figures/adr_', label, '1-*.png +append +repage -bordercolor White -border 3x3 figures/adr_', label, '1.png'];
    c2 = ['magick figures/adr_', label, '2-*.png +append +repage -bordercolor White -border 3x3 figures/adr_', label, '2.png'];
    c3 = ['magick figures/adr_', label, '1.png figures/adr_', label, '2.png -append figures/adr_', label, '.png'];
    c4 = ['rm figures/adr_', label, '[12]*.png'];
    system(c1);
    system(c2);
    system(c3);
    system(c4);
end
