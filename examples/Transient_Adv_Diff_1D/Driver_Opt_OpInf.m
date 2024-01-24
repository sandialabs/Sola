clear;
close all;
clc;
addpath(genpath('../../src'));

%% Experimental parameters.

datafile = 'OpInf_Training_Data.mat';
residual_energy_threshold = 1e-5;
plot_basis_functions = false;
plot_optimization_solution = false;

%% Generate or load training data.

if ~exist(datafile, 'file')
    % Set up the optimization problem.
    n_y = 200;
    n_t = 151;
    T = 1;
    num_space_control_nodes = 10;
    n_z = num_space_control_nodes * (n_t - 1);
    con = Adv_Diff_Gaussian_Source_Constraint(n_y, n_z, T, n_t, num_space_control_nodes);

    % Solve the state equation for several random controls.
    num_solves = 10;
    Z_train = randn(n_z, num_solves);
    U_train = zeros(n_y * n_t, num_solves);
    for k = 1:num_solves
        U_train(:, k) = con.State_Solve(Z_train(:, k));
    end

    % Save the experimental parameters and the training data.
    save(datafile, "n_y", "n_z", "T", "n_t", "num_space_control_nodes", "num_solves", "Z_train", "U_train");
else
    load(datafile);
end

%% Set up objectives and learn POD basis.

obj_hifi = Adv_Diff_Gaussian_Source_Objective(n_y, n_z, T, n_t, num_space_control_nodes);
states = reshape(U_train, n_y, n_t * num_solves);
basis = POD_Basis(states, false, obj_hifi.M);
basis.Set_Reduced_Dimension_From_Residual_Energy(residual_energy_threshold);
obj_lofi = Reduced_Dynamic_Objective(obj_hifi, basis.V);
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

% Check projection error.
states_lofi = basis.Compress(states);
states_projected = basis.Decompress(states_lofi);
state_projection_error = norm(states_projected - states) / norm(states);
disp(['Projection error of training states: ', num2str(state_projection_error)]);

target_projection_error = zeros(1, n_t - 1);
for k = 2:n_t
    target = obj_hifi.Evaluate_Target(obj_hifi.t_mesh(k), obj_hifi.x);
    target_projected = basis.Project(target);
    target_projection_error(k) = norm(target - target_projected) / norm(target);
end
disp(['Average projection error of target state: ', num2str(mean(target_projection_error))]);

%% Initialize and calibrate Operator Inference constraint.

operators = {Linear_Operator(), Input_Operator()};
rom = OpInf_ROM_Constraint(basis.r, num_space_control_nodes, T, n_t, zeros(basis.r, 1), operators);
Qhats = reshape(states_lofi, basis.r * n_t, num_solves);
rom.Select_Regularization(Qhats, Z_train);

% Validate the Operator Inference constraint by solving
% the ROM for each of the training controls.
for k = 1:num_solves
    u = U_train(:, k);
    z = Z_train(:, k);
    u_rom = reshape(basis.Decompress(reshape(rom.State_Solve(z), basis.r, n_t)), n_y * n_t, 1);
    state_error = norm(u - u_rom) / norm(u);
    disp(['ROM reconstruction error for training set ', num2str(k), ': ', num2str(state_error)]);
end

%% Solve the optimization problem.

opt = Reduced_Space_Optimization(obj_lofi, rom);
z0 = rand(n_z, 1);
[u_reduced, z] = opt.Optimize(z0);

%% Compare the optimal state to the target.

x = obj_hifi.x;
t = obj_lofi.t_mesh;

u_reshape = reshape(u_reduced, basis.r, n_t);

relative_errors = zeros(n_t, 1);
if plot_optimization_solution
    figure;
end
for k = 1:n_t
    target = obj_hifi.Evaluate_Target(t(k), x);
    current = basis.Decompress(u_reshape(:, k));
    relative_errors(k) = norm(target - current) / norm(target);
    if plot_optimization_solution
        plot(x, current, '-', x, target, '--', 'LineWidth', 3);
        legend({'State', 'Target'});
        ylim([0 .2]);
        pause(.05);
    end
end
