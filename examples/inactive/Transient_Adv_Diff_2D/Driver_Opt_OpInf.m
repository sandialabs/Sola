%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
close all;
clc;

%% Experiment parameters.

meshfile = 'urban_canyon.mat';
datafile = 'OpInf_Training_Data.mat';
regenerate_data = false;
residual_energy_threshold = 1e-2;
plot_basis_functions = true;
plot_training_data = false;
plot_training_reconstruction = false;
basis_sizes = [];  % 2:4:40;
regularization_candidates = logspace(-8, 3, 80);
ddt_strategy = '6thOrder';

%% Generate training data if needed.

if ~exist(datafile, 'file') || regenerate_data
    disp('Generating training data');

    tic();
    % Initial condition parameters.
    init_centers = [.05 .85
                    .10 .50
                    .60 .70
                    1.1 .80
                    .80 .20]';
    num_solves = size(init_centers, 2);

    % Input function parameters.
    control_nodes = [0.1 0.5
                     0.1 0.9
                     0.1 1.1
                     0.3 0.7
                     0.3 0.9
                     0.3 1.1
                     0.5 0.3
                     0.5 0.5
                     0.5 0.7
                     0.7 0.7
                     0.9 0.3
                     0.9 1.1
                     1.1 0.7
                     1.1 0.9]';
    n_q = size(control_nodes, 2);

    % Time domain.
    t = linspace(0, .4, 101);
    n_t = length(t);
    n_z = (n_t - 1) * n_q;

    % Load spatial geometry and mesh.
    model = Transient_Adv_Diff_2D.model_fromfile(meshfile);
    n_y = size(model.Mesh.Nodes, 2);
    n_u = n_y * n_t;

    % Model and input parameters.
    diffusion = 0.005;
    advection = 4;
    num_randcontrol_nodes = 3;
    randcontrol_nodes = linspace(t(1), t(end), num_randcontrol_nodes);

    Z_train = zeros(n_z, num_solves);
    U_train = zeros(n_u, num_solves);

    for k = 1:num_solves
        disp(['High-fidelity solve ', num2str(k)]);

        % Initialize the solver.
        solver = Transient_Adv_Diff_2D(model, init_centers(:, k), ...
                                       diffusion, advection, control_nodes);

        % Set up a random control profile.
        vals = [zeros(n_q, 1), -20 * rand(n_q, num_randcontrol_nodes - 1)];
        pp = spline(randcontrol_nodes, vals);
        controller = @(tt) ppval(pp, tt);

        % Solve the system.
        Yk = solver.State_Solve(controller, t, plot_training_data).NodalSolution;
        Qk = controller(t(2:end));

        % Record results.
        U_train(:, k) = reshape(Yk, [], 1);
        Z_train(:, k) = reshape(Qk, [], 1);
    end
    time_trainingdata = toc();

    save(datafile, "t", "solver", "U_train", "Z_train", "time_trainingdata");
end

%% Load training data.

load(datafile);
n_t = length(t);
T = t(end);
n_u = size(U_train, 1);
num_solves = size(U_train, 2);
n_y = n_u / n_t;  % = size(solver.model.Mesh.Nodes, 2);
mass_matrix = assembleFEMatrices(solver.model, 'M').M;
n_z = size(Z_train, 1);
n_q = n_z / (n_t - 1);  % = solver.n_q;
disp(['Using ', num2str(num_solves), ' training trajectories']);

%% Learn a POD basis.

% Unpack the states and controls by training trajectory.
states = cell(num_solves);
controls = cell(num_solves);
for k = 1:num_solves
    states{k} = reshape(U_train(:, k), n_y, n_t);
    controls{k} = reshape(Z_train(:, k), n_q, n_t - 1);
end

% Learn a POD basis from the collection of all state snapshots.
states_all = horzcat(states{:});
basis = POD_Basis(states_all, false, full(mass_matrix));
basis.Set_Reduced_Dimension_From_Residual_Energy(residual_energy_threshold);
disp(['Selected reduced dimension r = ', num2str(basis.r)]);

if plot_basis_functions
    for j = 1:basis.r
        solver.Plot_Field(basis.V(:, j));
        title(['POD basis function ', num2str(j)]);
    end
end

if plot_training_data
    for k = 1:num_solves
        Yhatk = basis.Compress(states{k});
        figure;
        plot(t, Yhatk);
        title(['compressed state training data, trajectory', num2str(k)]);
    end
end

%% Learn a ROM, varying the reduced state dimension.

if isempty(basis_sizes)
    basis_sizes = [basis.r];
end
errors = zeros(length(basis_sizes), 1);
for i = 1:length(basis_sizes)
    r = basis_sizes(i);
    basis.r = r;
    disp(' ');
    disp(['POD with ' num2str(r), ' basis vectors']);

    % Compress states and check projection error.
    states_lofi = cell(num_solves);
    for k = 1:num_solves
        states_lofi{k} = basis.Compress(states{k});
        Yk_proj = basis.Decompress(states_lofi{k});
        proj_err_k = norm(Yk_proj - states{k}) / norm(states{k});
        disp(['Projection error of training states for trajectory ', ...
              num2str(k), ': ', ...
              num2str(100 * proj_err_k), '%']);
    end

    %% Learn an OpInf ROM from the data.

    operators = {Linear_Operator(), Input_Operator()};
    rom = OpInf_ROM_Constraint(basis.r, n_q, T, n_t, zeros(basis.r, 1), operators);
    tic();
    rom.Select_Regularization(states_lofi, controls, regularization_candidates, ddt_strategy);
    time_opinfcalibration = toc();

    % Solve the ROM for each of the training controls.
    total_error = 0;
    for k = 1:num_solves
        Yk_data = states{k};
        rom.y0 = states_lofi{k}(:, 1);
        Yk_rom = basis.Decompress(rom.State_Solve2(controls{k}));
        state_error = norm(Yk_data - Yk_rom) / norm(Yk_data);
        disp(['ROM reconstruction error for training set ', num2str(k), ': ', num2str(100 * state_error), '%']);
        total_error = total_error + state_error;
        if plot_training_reconstruction
            solver.Animate_Solution(Yk_rom);
        end
    end
    errors(i) = total_error / num_solves;
end

if length(basis_sizes) > 1
    figure;
    plot(basis_sizes, errors);
    title('POD basis vector size versus average ROM training error');
end

%% Set up and solve the optimization problem.

% Make sure the initial conditions are right.
solver.init_center = [.05; .85];
rom.y0 = states_lofi{1}(:, 1);

obj_hifi = solver.Make_Objective([.6; .6], t(end), length(t), 1e-5);
obj_lofi = Reduced_Dynamic_Objective(obj_hifi, basis.V);
solver.Plot_Field(obj_hifi.target_weight, 'Protection zone');

opt = Reduced_Space_Optimization(obj_lofi, rom);
opt.z_lb = -100 * ones(n_z, 1);             % Lower bounds for control.
opt.z_ub = zeros(n_z, 1);                   % Upper bounds for control.

tic();
[u_lofi, z_lofi] = opt.Optimize(-rand(n_z, 1));
time_lofioptimization = toc();

% Inspect the state solution.
Y_rom = basis.Decompress(reshape(u_lofi, basis.r, n_t));
solver.Animate_Solution(Y_rom);             % ROM state with ROM controller

% Inspect the control solution.
Q_rom = reshape(z_lofi, n_q, n_t - 1);
figure;
plot(t(2:end), Q_rom);
title('Optimal controls (optimized with an OpInf ROM surrogate)');

% Solve the high-fidelity model with the inferred controls.
disp('Final high-fidelity solve');
pp = spline(t(2:end), Q_rom);
controller = @(tt) ppval(pp, tt);
Y_hifi = solver.State_Solve(controller, t, false).NodalSolution;
solver.Animate_Solution(Y_hifi);            % FOM state with ROM controller

save('OptimizationSolution.mat', "solver", "Y_hifi", "Y_rom", "t", "Q_rom", "n_q");

%% Load and visualize results later.
% load('OptimizationSolution.mat', "solver", "Y_hifi", "t", "z_lofi", "n_q");
% figure;
% plot(t(2:end), Q_rom);
% title('Optimal controls (optimized with an OpInf ROM surrogate)');
% solver.Animate_Solution(Y_rom);   % ROM state with ROM controller
% solver.Animate_Solution(Y_hifi);  % FOM state with ROM controller
