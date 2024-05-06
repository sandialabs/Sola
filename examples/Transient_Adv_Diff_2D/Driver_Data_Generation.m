clear;
close all;
clc;
run('../../src/Set_Paths');

% Beforehand: use pdeModeler to generate geometry and mesh.
% save('urban_canyon.mat', 'points', 'edges', 'triangles');

control_nodes = [
                 % 0.1000    0.1000
                 % 0.1000    0.3000
                 0.1000    0.5000
                 % 0.1000    0.7000
                 0.1000    0.9000
                 0.1000    1.1000
                 % 0.3000    0.1000
                 % 0.3000    0.3000
                 % 0.3000    0.5000
                 0.3000    0.7000
                 0.3000    0.9000
                 0.3000    1.1000
                 % 0.5000    0.1000
                 0.5000    0.3000
                 0.5000    0.5000
                 0.5000    0.7000
                 % 0.5000    0.9000
                 % 0.5000    1.1000
                 % 0.7000    0.1000
                 % 0.7000    0.3000
                 % 0.7000    0.5000
                 0.7000    0.7000
                 % 0.7000    0.9000
                 % 0.7000    1.1000
                 % 0.9000    0.1000
                 0.9000    0.3000
                 % 0.9000    0.5000
                 % 0.9000    0.7000
                 % 0.9000    0.9000
                 0.9000    1.1000
                 % 1.1000    0.1000
                 % 1.1000    0.3000
                 % 1.1000    0.5000
                 1.1000    0.7000
                 1.1000    0.9000
                 % 1.1000    1.1000
                ]';

%% Initialize the solver.
% model = Transient_Adv_Diff_2D.model_default(.05);
model = Transient_Adv_Diff_2D.model_fromfile('urban_canyon.mat');
solver = Transient_Adv_Diff_2D(model, [.05; .85], .01, 4, control_nodes);
t = linspace(0, .2, 400);

%% Visualize the solver geometry.
solver.Plot_Control_Nodes();
solver.Plot_Velocity_Field();

%% Initialize the controller.
controller = randomspline(t, solver.n_q, 4);
% controller = sinecontrol(solver.n_q);
% controller = nocontrol(solver.n_q);

%% Solve with the selected controller.
u = solver.State_Solve(controller, t, true);

%% Plot the initial condition.
solver.Plot_Field(u.NodalSolution(:, 1), 'Initial condition');

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
    vals = -10 * rand(dofs, num_nodes);
    pp = spline(nodes, vals);
    out = @(tt) ppval(pp, tt);
end

function [out] = sinecontrol(dofs)
    % Sinusoidal inputs.
    out = @(t) 10 * sin(16 * pi * t) .* ones(dofs, 1);
end

function [out] = nocontrol(dofs)
    % Zero inputs: q(t) = 0 for all t/
    out = @(t) zeros(dofs, length(t));
end
