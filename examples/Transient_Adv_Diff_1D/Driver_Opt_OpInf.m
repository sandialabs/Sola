clear;
close all;
clc;
addpath(genpath('../../src'));

%% Experimental parameters.
datafile = 'OpInf_Training_Data.mat';
residual_energy_threshold = 1e-5;
plot_basis_functions = false;
plot_optimization_solution = false;

%% Load data.
if ~exist(datafile, 'file')
    run Driver_Data_Generation.m;
end
load(datafile);

%% Set up objectives and learn POD basis.
obj_hifi = Adv_Diff_Gaussian_Source_Objective(n_y, n_z, T, n_t, num_space_control_nodes);
states = reshape(Y, n_y, n_t * num_solves);
basis = POD_Basis(states, false, obj_hifi.M);
basis.Set_Reduced_Dimension_From_Residual_Energy(residual_energy_threshold);
disp(['Selected reduced dimension r = ', num2str(basis.r)]);

if plot_basis_functions
    figure;
    plot(obj_hifi.x, basis.V(:, 1));
    hold on;
    for j = 2:basis.r
        plot(obj_hifi.x, basis.V(:, j));
    end
    title('POD basis functions');
end

obj_lofi = Reduced_Dynamic_Objective(obj_hifi, basis.V);

%% Compress states and estimate time derivatives of states.
states_lofi = basis.Compress(states);

Qhats = reshape(states_lofi, basis.r, n_t, num_solves);
ddts = zeros(basis.r, n_t - 1, num_solves);
dt = obj_hifi.t_mesh(2) - obj_hifi.t_mesh(1);
for k = 1:num_solves
    ddts(:, :, k) = Qhats(:, 2:end, k) - Qhats(:, 1:(end - 1), k) / dt;
end
states_lofi = reshape(Qhats(:, 1:(end - 1), :), basis.r, (n_t - 1) * num_solves);
ddts = reshape(ddts, basis.r, (n_t - 1) * num_solves);
inputs = reshape(Z, num_space_control_nodes, (n_t - 1) * num_solves);
y0 = zeros(basis.r, 1);

%% Initialize Operator Inference constraint.

operators = {Linear_Operator(), Input_Operator()};
rom = OpInf_ROM_Constraint(basis.r, num_space_control_nodes, T, n_t, y0, operators);
rom.Learn_Operators(states_lofi, inputs, ddts);

error('Done to here');

%% Solve the optimization problem.
opt = Reduced_Space_Optimization(obj_lofi, rom);
z0 = rand(n_z, 1);
[u, z] = opt.Optimize(z0);

%% Compare the optimal state to the target.
x = obj_hifi.x;
t = obj_lofi.t_mesh;

u_reshape = reshape(u, n_y, n_t);

if plot_optimization_solution
    figure;
    for k = 1:n_t
        target = obj_hifi.Evaluate_Target(t(k), x);
        plot(x, u_reshape(:, k), '-', x, target, '--', 'LineWidth', 3);
        legend({'State', 'Target'});
        ylim([0 .2]);
        pause(.05);
    end
end

% TODO: error analysis.
