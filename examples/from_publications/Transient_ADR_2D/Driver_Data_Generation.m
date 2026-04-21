%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
close all;
clc;

% Beforehand: use pdeModeler to generate geometry and mesh.
% save('urban_canyon.mat', 'points', 'edges', 'triangles');

control_nodes = [
                 0.1000    0.5000
                 0.1000    0.9000
                 0.1000    1.1000
                 0.3000    0.7000
                 0.3000    0.9000
                 0.3000    1.1000
                 0.5000    0.3000
                 0.5000    0.5000
                 0.5000    0.7000
                 0.7000    0.7000
                 0.9000    0.3000
                 0.9000    1.1000
                 1.1000    0.7000
                 1.1000    0.9000
                ]';

%% Initialize the solver.
% model = Transient_ADR_2D.model_default(.02);
% solver = Transient_ADR_2D(model, [-0.5; 0.5], [.05; .05], [4; 4], 10, 4);

model = Transient_ADR_2D.model_fromfile('urban_canyon.mat');
solver = Transient_ADR_2D(model, [.05; .85], [.1; .1], [4; 4], 2, control_nodes);

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

%% Solve with the selected controller.
tic();
u = solver.State_Solve(controller, t);
solve_time = toc();

%% Visualize the results.
solver.Plot_Field(u.NodalSolution(:, :, 1), 'Initial condition');
solver.Animate_Solution(u.NodalSolution);

%% Save the results.
save('solver.mat', 'solver');
save('solution.mat', 'u');

% % For later: load and animate the results again.
% load('solver.mat', 'solver');
% load('solution.mat', 'u');
% solver.Animate_Solution(u.NodalSolution)

%% Controllers
function [out] = randomspline(t, dofs, num_nodes)
    % Random (but smooth) nonpositive inputs.
    nodes = linspace(min(t), max(t), num_nodes);
    vals = 50 * rand(dofs, num_nodes) + 5;
    pp = spline(nodes, vals);
    out = @(tt) ppval(pp, tt);
end

function [out] = nocontrol(dofs)
    % Zero inputs: q(t) = 0 for all t/
    out = @(t) zeros(dofs, length(t));
end
