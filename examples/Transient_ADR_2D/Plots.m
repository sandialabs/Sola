clear;
close all;
clc;
run('../../src/Set_Paths');

% meshfile = 'urban_canyon.mat';
datafile = 'OpInf_Training_Data.mat';

%% Plot configuration.

set(0, 'DefaultAxesFontName', 'Times New Roman');
set(0, 'DefaultAxesFontWeight', 'normal');
set(0, 'DefaultAxesFontSize', 14);
set(0, 'DefaultTextFontName', 'Times New Roman');
set(0, 'DefaultTextFontWeight', 'normal');
set(0, 'DefaultTextFontSize', 18);

%% Load the model mesh and model.

% model = Transient_ADR_2D.model_fromfile(meshfile);
load(datafile, 't', 'solver');
load('fem_matrices.mat', 'mass_matrix');
objective = Transient_ADR_2D_Objective([.6; .6], solver.x, solver.y, mass_matrix, t(end), length(t), solver.n_q, 1e-4);

%% Finite element mesh, protection zone over the spatial domain.

% Finite element mesh.
fig1 = figure;
ax1 = subplot(1, 1, 1);
pdemesh(solver.model);
hold on;
scatter(ax1, solver.control_nodes(1, :), solver.control_nodes(2, :), 36, "black", "filled", "o");
xlim(ax1, [min(solver.x, [], "all") max(solver.x, [], "all")]);
ylim(ax1, [min(solver.y, [], "all") max(solver.y, [], "all")]);
t1 = title(ax1, "Finite Element Mesh and Injection Locations");
set(t1, 'FontWeight', 'normal');

% Protection zone.
fig2 = figure;
ax2 = subplot(1, 1, 1);
pdeplot(solver.model.Mesh, XYData = objective.target_weight, ColorMap = "parula");
c2 = colorbar(ax2);
xlim(ax2, [min(solver.x, [], "all") max(solver.x, [], "all")]);
ylim(ax2, [min(solver.y, [], "all") max(solver.y, [], "all")]);
t2 = title(ax2, "Protection weights \psi");
set(t2, 'FontWeight', 'normal');

% Format axes and figure.
linkaxes([ax1, ax2], 'xy');
axis(ax2, 'equal');
for ax = [ax1 ax2]
    box(ax, 'off');
    set(ax, 'XTick', [], 'YTick', [], 'XTickLabel', [], 'YTickLabel', []);
    set(ax, 'Color', 'none');
    set(ax, 'XColor', 'none', 'YColor', 'none');
end
pos1 = get(ax1, 'Position');
pos2 = get(ax2, 'Position');
pos2(3:4) = pos1(3:4);
set(ax2, 'Position', pos2);
for fig = [fig1 fig2]
    set(fig, 'PaperPositionMode', 'auto');
    set(fig, 'Units', 'inches');
    set(fig, 'PaperUnits', 'inches');
    figPosition = get(fig, 'Position');
    figWidth = figPosition(3);
    figHeight = figPosition(4);
    set(fig, 'PaperSize', [figWidth figHeight]);
    set(fig, 'PaperPosition', [0 0 figWidth figHeight]);
end

% Save the figures.
print(fig1, 'figures/adr_mesh.pdf', '-dpdf', '-r300', '-loose');
print(fig2, 'figures/adr_target.png', '-dpng', '-r300', '-loose');
% Postprocessing (on MacOS):
%  - Download TeX from https://www.tug.org/mactex/mactex-download.html
%  - Install brew from https://brew.sh/
%  - Run the following commands in a terminal.
%    $ pdfcrop --margin "1" figures/adr_mesh.pdf figures/adr_mesh.pdf
%    $ brew install imagemagick     # Only need to do this once ever.
%    $ magick adr_target.png -trim adr_target.png

%% POD singular value decay for each state variable.

% TODO

%% ROMCO control solution as a function of time.

% TODO

%% FOM solution with ROMCO controls as a function of time (including init).

% TODO
