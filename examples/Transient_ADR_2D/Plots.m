clear;
close all;
clc;
run('../../src/Set_Paths');

% meshfile = 'urban_canyon.mat';
datafile = 'OpInf_Training_Data.mat';

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
load(datafile, 't', 'solver');
load('fem_matrices.mat', 'mass_matrix');
objective = Transient_ADR_2D_Objective([.6; .6], solver.x, solver.y, mass_matrix, t(end), length(t), solver.n_q, 1e-4);
load('bases.mat', 'svdvals', 'svdvecs1', 'svdvecs2');

%% Finite element mesh.
fig = plotfield(solver, "mesh", 0, true, "Spatial mesh and injection locations");
print(fig, 'figures/adr_mesh.pdf', '-dpdf', '-r300', '-loose');
close(fig);

%% Protection zone.
fig = plotfield(solver, objective.target_weight, viridis, true, "Protection zone weights \psi");
print(fig, 'figures/adr_target.png', '-dpng', '-r300', '-loose');
close(fig);

%% Initial condition.
u10 = 50 * exp(-50 .* sum(([solver.x, solver.y]' - solver.init_center).^2, 1))';
fig = plotfield(solver, u10, viridis, true, "Initial condition for contaminant");
print(fig, 'figures/adr_initial.png', '-dpng', '-r300', '-loose');
close(fig);

%% POD singular value decay for each state variable.
% j = 0:size(svdvals, 2);
% resenergy = 1 - (cumsum(svdvals.^2, 2) ./ sum(svdvals.^2, 2));
% resenergy = [1, resenergy(1, :); 1, resenergy(2, :)];
%
% fig = figure;
% ax = subplot(1, 1, 1);
% for i = [1, 2]
%     semilogy(ax, j, resenergy(i, :), '.-', MarkerSize=10);
%     hold on
% end
% xline(ax, 24, 'Color', ax4.ColorOrder(1, :));
% xline(ax, 30, 'Color', ax4.ColorOrder(2, :));
% yline(ax, 1e-5);
% xlim(ax, [0, 50]);
% ylim(ax, [1e-8, 1e0]);
% xlabel(ax, 'Singular value index');
% ylabel(ax, 'Residual energy');
% text(ax, 23.5, 2e-8, '$r_1 = 24$', 'Interpreter', 'latex', ...
%      'HorizontalAlignment', 'right', 'Color', ax4.ColorOrder(1, :), ...
%      'FontSize', 16);
% text(ax, 30.5, 2e-8, '$r_2 = 30$', 'Interpreter', 'latex', ...
%      'HorizontalAlignment', 'left', 'Color', ax4.ColorOrder(2, :), ...
%      'FontSize', 16);
% print(fig, 'figures/adr_svdvals.pdf', '-dpdf', '-r300', '-loose');

%% Dominant basis vectors for each variable.
% [fig, ax] = plotfield(solver, svdvecs1(:, 1), viridis, true, "First basis function for contaminant");
% print(fig, 'figures/adr_basis1.png', '-dpng', '-r300', '-loose');
% close(fig);
%
% [fig, ax] = plotfield(solver, svdvecs2(:, 1), "parula", true, "First basis function for neutralizer");
% print(fig, 'figures/adr_basis2.png', '-dpng', '-r300', '-loose');
% close(fig);

%% Grid of basis vectors over the spatial domain.
lim1 = max(abs([min(svdvecs1(:)), max(svdvecs1(:))]));
lim2 = max(abs([min(svdvecs2(:)), max(svdvecs2(:))]));
for j = 1:4
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

% TODO

%% FOM solution with ROMCO controls as a function of time (including init).

% TODO

%% Postprocess saved figures.

[os, ~] = computer;
if strcmp(os, 'MACA64') || strcmp(os, 'MACI64')
    % Postprocessing (on MacOS):
    %  - Download TeX from https://www.tug.org/mactex/mactex-download.html
    %  - Install brew from https://brew.sh/
    %  - Check `which pdfcrop` and `which magick` give the paths below.
    setenv('PATH', [getenv('PATH') ':/Library/TeX/texbin' ':/opt/homebrew/bin']);
    !for fname in figures/*.pdf; do pdfcrop --margin "1" "${fname}" "${fname}"; done
    !for fname in figures/*.png; do magick "${fname}" -trim +repage -bordercolor White -border 10x20 "${fname}"; done
    !magick figures/adr_basis1-*.png +append +repage -bordercolor White -border 3x3 figures/adr_basis1.png
    !magick figures/adr_basis2-*.png +append +repage -bordercolor White -border 3x3 figures/adr_basis2.png
    !magick figures/adr_basis1.png figures/adr_basis2.png -append figures/adr_basis.png
    !rm figures/adr_basis[12]*.png
    disp('ALL FIGURES PROCESSED AND SAVED');
else
    disp('Postprocessing only defined for MacOS');
end

%% Helper functions

function [fig, ax] = plotfield(solver, data, cmp, cb, titletext, ylabeltext, limits)
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
    end

    fig = figure;
    ax = subplot(1, 1, 1);

    % Plot the data, add colorbar, etc.
    if isstring(data) && data == "mesh"
        pdemesh(solver.model);
        hold on;
        scatter(ax, solver.control_nodes(1, :), solver.control_nodes(2, :), 36, "black", "filled", "o");
    else
        pdeplot(solver.model.Mesh, XYData = data, ColorMap = cmp, ColorBar = "off");
        if cb
            % Colorbar.
            cbar = colorbar(ax);
            colorbarPos = cbar.Position;
            colorbarPos(1) = colorbarPos(1) - 0.01;
            cbar.Position = colorbarPos;
        end
        if ~isempty(limits)
            clim(ax, limits);
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
