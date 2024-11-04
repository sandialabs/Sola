clear;
close all;
clc;
run('../../src/Set_Paths');

% Beforehand: use pdeModeler to generate geometry and mesh.
% save('mesh.mat', 'geometry', 'bcs', 'points', 'edges', 'triangles');

%% Define control nodes.
control_grid_size = 9;
control_locs = linspace(0, 1.2, control_grid_size + 2);
control_locs = control_locs(2:end - 1);
control_nodes = table2array(combinations(control_locs, control_locs))';

%% Initialize the solver.
model = Transient_ADR_2D.model_fromfile('mesh.mat');
solver = Transient_ADR_2D(model, control_nodes);

t = linspace(0, .4, 101);

%% Visualize the solver geometry.
solver.Plot_Control_Nodes();
solver.Plot_Velocity_Field();

%% Initialize the controller and visualize the controls.
controller = randomspline(t, solver.n_q, 8);
% controller = nocontrol(solver.n_q);

figure;
plot(t, controller(t));
title('Controller');

% solver.Animate_Controls(controller(t));

%% Solve the model with the selected controller.
tic();
u = solver.State_Solve(controller, t);
solve_time = toc();

%% Visualize the results.
solver.Plot_Field(u.NodalSolution(:, :, 1), 'Initial condition');
solver.Animate_Solution(u.NodalSolution, false);

%% Contour animation.
solver.Animate_Contours(u.NodalSolution);

%% Save the results.
save('fom_solver.mat', 'solver');
save('fom_solution.mat', 'u');

% % For later: load and animate the results again.
% load('fom_solver.mat', 'solver');
% load('fom_solution.mat', 'u');
% solver.Animate_Solution(u.NodalSolution);

%% Controllers
function [out] = randomspline(t, dofs, num_nodes)
    % Random (but smooth) nonnegative inputs.
    nodes = linspace(min(t), max(t), num_nodes);
    vals = 20 * rand(dofs, num_nodes);
    pp = pchip(nodes, vals);
    out = @(tt) ppval(pp, tt);
end

function [out] = nocontrol(dofs)
    % Zero inputs: q(t) = 0 for all t/
    out = @(t) zeros(dofs, length(t));
end
