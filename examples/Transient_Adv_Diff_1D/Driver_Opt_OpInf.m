clear;
close all;
clc;
addpath(genpath('../../src'));

%% Experimental parameters.

datafile = 'OpInf_Training_Data.mat';
residual_energy_threshold = 1e-5;
plot_basis_functions = true;
plot_optimization_solution = true;

%% Generate or load training data.

if ~exist(datafile, 'file')
    tic();
    % Set up the optimization problem.
    n_y = 200;
    n_t = 51;
    T = 1;
    num_space_control_nodes = 10;
    n_z = num_space_control_nodes * (n_t - 1);
    con = Adv_Diff_Gaussian_Source_Constraint(n_y, n_z, T, n_t, num_space_control_nodes);

    % Solve the state equation for several random controls.
    num_solves = 3;
    Z_train = randn(n_z, num_solves);
    U_train = zeros(n_y * n_t, num_solves);
    for k = 1:num_solves
        U_train(:, k) = con.State_Solve(Z_train(:, k));
    end
    time_trainingdata = toc();

    % Save the experimental parameters and the training data.
    save(datafile, "n_y", "n_z", "T", "n_t", "num_space_control_nodes", "num_solves", "Z_train", "U_train", "time_trainingdata");
else
    load(datafile);
end
disp(['Using ', num2str(num_solves), ' training trajectories']);

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
tic();
rom.Select_Regularization(Qhats, Z_train, logspace(-8, -3, 20));
time_opinfcalibration = toc();

% Validate the Operator Inference constraint by solving
% the ROM for each of the training controls.
for k = 1:num_solves
    u = reshape(U_train(:, k), n_y, n_t);
    z = Z_train(:, k);
    u_rom_k = basis.Decompress(reshape(rom.State_Solve(z), basis.r, n_t));
    state_error = norm(u - u_rom_k) / norm(u);
    disp(['ROM reconstruction error for training set ', num2str(k), ': ', num2str(state_error)]);
end

%% Solve the optimization problem.

opt = Reduced_Space_Optimization(obj_lofi, rom);
z0 = rand(n_z, 1);

tic();
[u_reduced, z_lofi] = opt.Optimize(z0);
time_lofioptimization = toc();

u_rom = basis.Decompress(reshape(u_reduced, basis.r, n_t));

%% For comparison, solve the optimization problem in full fidelity.

con_hifi = Adv_Diff_Gaussian_Source_Constraint(n_y, n_z, T, n_t, num_space_control_nodes);
opt_hifi = Reduced_Space_Optimization(obj_hifi, con_hifi);

tic();
[u_true, z_hifi] = opt_hifi.Optimize(z0);
time_hifioptimization = toc();

u_true = reshape(u_true, n_y, n_t);
diff_state = norm(u_rom - u_true) / norm(u_true);
diff_control = norm(z_hifi - z_lofi) / norm(z_hifi);

%% Compare the optimal state to the target.

x = obj_hifi.x;
t = obj_hifi.t_mesh;
err1 = zeros(n_t, 1);
err2 = zeros(n_t, 1);
if plot_optimization_solution
    figure;
end
for k = 1:n_t
    target = obj_hifi.Evaluate_Target(t(k), x);
    denom = norm(target);
    err1(k) = norm(target - u_rom(:, k)) / denom;
    err2(k) = norm(target - u_true(:, k)) / denom;
    if plot_optimization_solution
        plot(x, u_rom(:, k), '-', x, u_true(:, k), ':', x, target, '--', 'LineWidth', 3);
        legend({'State (ROM)', 'State (FOM)', 'Target'});
        ylim([0 .2]);
        pause(.05);
    end
end
if plot_optimization_solution
    figure;
    semilogy(t, err1, '-', t, err2, ':', 'LineWidth', 3);
    legend({'Target error (ROM)', 'Target error (FOM)'});
end

%% Final report.

disp('Error Report');
disp(['  Relative error of low-fi and high-fi opt state:', 9, num2str(diff_state)]);
disp(['  Relative error of low-fi and high-fi opt control:' 9, num2str(diff_control)]);
disp('Timing Report');
disp(['  Generate ', num2str(num_solves), ' training trajectories:', 9, num2str(time_trainingdata), 's']);
disp(['  OpInf ROM calibration:', 9, 9, num2str(time_opinfcalibration), 's']);
disp(['  OpInf-constrained optimization:', 9, num2str(time_lofioptimization), 's']);
disp(['  FOM-constrained optimization:', 9, 9, num2str(time_hifioptimization), 's']);
