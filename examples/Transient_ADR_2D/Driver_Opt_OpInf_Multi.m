clear;
close all;
clc;
run('../../src/Set_Paths');

%% Experiment parameters.

meshfile = 'urban_canyon.mat';
datafile = 'OpInf_Training_Data.mat';

regenerate_data = false;
plot_basis_functions = false;
plot_training_data = false;
plot_training_reconstruction = false;

residual_energies = [1e-5];
ABregularization_candidates = [1e-6, 1e-2, 1e1];
Hregularization_candidates = logspace(3, 5, 21);
ddt_strategy = '6thOrder';
control_regularization = 1.e-4;

%% Generate training data if needed.

if ~exist(datafile, 'file') || regenerate_data
    disp('Generating training data');

    tic();
    % Initial condition parameters.
    init_center = [.05; .85];
    num_solves = 5;

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
    model = Transient_ADR_2D.model_fromfile(meshfile);
    n_x = size(model.Mesh.Nodes, 2);
    n_y = 2 * n_x;
    n_u = n_y * n_t;

    % Model and input parameters.
    diffusion = [0.10, 0.10];
    advection = [4.00, 4.00];
    reaction = 2;
    num_randcontrol_nodes = 4;
    randcontrol_nodes = linspace(t(1), t(end), num_randcontrol_nodes);

    Z_train = zeros(n_z, num_solves);
    U_train = zeros(n_u, num_solves);

    for k = 1:num_solves
        disp(['High-fidelity solve ', num2str(k)]);

        % Initialize the solver.
        solver = Transient_ADR_2D(model, init_center, ...
                                  diffusion, advection, reaction, control_nodes);

        % Set up a random control profile.
        vals = [zeros(n_q, 1), 50 * rand(n_q, num_randcontrol_nodes - 1)];
        pp = spline(randcontrol_nodes, vals);
        controller = @(tt) ppval(pp, tt);

        % Solve the system.
        Yk = solver.State_Solve(controller, t).NodalSolution;

        if plot_training_data
            solver.Animate_Solution(Yk);
        end

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
n_y = n_u / n_t;
n_x = n_y / 2;  % = size(solver.model.Mesh.Nodes, 2);
mass_matrix = assembleFEMatrices(solver.model, 'M').M;
mass_matrix = mass_matrix(1:n_x, 1:n_x);
stiffness_matrix = assembleFEMatrices(solver.model, 'K').K;
stiffness_matrix = stiffness_matrix(1:n_x, 1:n_x);

n_z = size(Z_train, 1);
n_q = n_z / (n_t - 1);  % = solver.n_q;
fprintf('Using %d training trajectories\n', num_solves);

%% Learn a POD basis for each variable.

% Unpack the states and controls by training trajectory.
states = cell(num_solves);
controls = cell(num_solves);
controls_romtraining = cell(num_solves);
for k = 1:num_solves
    states{k} = reshape(U_train(:, k), n_y, n_t);
    controls{k} = reshape(Z_train(:, k), n_q, n_t - 1);
    controls_romtraining{k} = sqrt(abs(controls{k}));
end

% Learn POD bases from the collection of all state snapshots.
states_all = horzcat(states{:});
basis1 = POD_Basis(states_all(1:n_x, :), false);  % , full(mass_matrix));
basis1.Set_Reduced_Dimension_From_Residual_Energy(residual_energies(1));
basis2 = POD_Basis(states_all(n_x + 1:end, :), false);  % , full(mass_matrix));
basis2.Set_Reduced_Dimension_From_Residual_Energy(residual_energies(1));

if plot_basis_functions
    for i = 1:min(basis1.r, basis2.r)
        solver.Plot_Field([basis1.V(:, i), basis2.V(:, i)]);
        title(['POD basis function ', num2str(i)]);
    end
end

if plot_training_data
    for k = 1:num_solves
        Yhatk_1 = basis1.Compress(states{k}(1:n_x, :));
        Yhatk_2 = basis2.Compress(states{k}(n_x + 1:end, :));
        Yhatk = [Yhatk_1; Yhatk_2];
        figure;
        plot(t, Yhatk);
        title(['compressed state training data, trajectory', num2str(k)]);
    end
end

%% Learn a ROM, varying the reduced state dimension.

errors = zeros(length(residual_energies), 1);
for i = 1:length(residual_energies)
    res_energy = residual_energies(i);
    fprintf('\nUsing %.2e residual energy\n', res_energy);

    basis1.Set_Reduced_Dimension_From_Residual_Energy(res_energy);
    basis2.Set_Reduced_Dimension_From_Residual_Energy(res_energy);
    r_1 = basis1.r;
    r_2 = basis2.r;
    n_yr = r_1 + r_2;
    fprintf('POD with r_1 = %d and r_2 = %d basis vectors\n', r_1, r_2);

    % Compress states and check projection error.
    states_lofi = cell(num_solves);
    for k = 1:num_solves
        Yhat_1 = basis1.Compress(states{k}(1:n_x, :));
        Yhat_2 = basis2.Compress(states{k}(n_x + 1:end, :));
        states_lofi{k} = [Yhat_1; Yhat_2];
        Yproj_1 = basis1.Decompress(Yhat_1);
        Yproj_2 = basis2.Decompress(Yhat_2);
        Yproj = [Yproj_1; Yproj_2];
        proj_err = norm(Yproj - states{k}) / norm(states{k});
        fprintf("Projection error of trajectory %d: %.4f%%\n", k, 100 * proj_err);
    end

    % Learn an OpInf ROM from the data.
    rom = Transient_ADR_2D_OpInf_Constraint(r_1, r_2, n_q, T, n_t, zeros(n_yr, 1));
    tic();
    rom.Select_Regularization(states_lofi, controls_romtraining, ...
                              ABregularization_candidates, ...
                              Hregularization_candidates, ...
                              ddt_strategy);
    time_opinfcalibration = toc();

    % Solve the ROM for each of the training controls.
    total_error = 0;
    for k = 1:num_solves
        Yk_data = states{k};
        rom.y0 = states_lofi{k}(:, 1);
        Yk_rom_compressed = rom.State_Solve2(controls_romtraining{k});
        Yk_rom_1 = basis1.Decompress(Yk_rom_compressed(1:r_1, :));
        Yk_rom_2 = basis2.Decompress(Yk_rom_compressed(r_1 + 1:end, :));
        Yk_rom = [Yk_rom_1; Yk_rom_2];
        state_error = norm(Yk_data - Yk_rom) / norm(Yk_data);
        fprintf('ROM reconstruction error for training set %d: %.2f%%\n', k, 100 * state_error);
        total_error = total_error + state_error;
        if plot_training_reconstruction
            solver.Animate_Solution(Yk_rom);
        end
    end
    errors(i) = total_error / num_solves;
end

if length(residual_energies) > 1
    figure;
    semilogx(residual_energies, errors);
    title('Residual energy versus average ROM training error');
end

%% Set up the optimization objective.

% Make sure the initial conditions are right.
solver.init_center = [.05; .85];
rom.y0 = states_lofi{1}(:, 1);

obj_hifi = solver.Make_Objective([.6; .6], t(end), length(t), control_regularization);
Vfull = blkdiag(basis1.V, basis2.V);
obj_lofi = Reduced_Dynamic_Objective(obj_hifi, Vfull);
solver.Plot_Field(obj_hifi.target_weight, 'Protection zone');

%% Set up and solve the optimization problem (using last trained ROM).

opt = Reduced_Space_Optimization(obj_lofi, rom);
% opt.z_lb = zeros(n_z, 1);                   % Lower bounds for control.
% opt.z_ub = 25 * ones(n_z, 1);               % Upper bounds for control.
opt.max_cg_iter = 200;

tic();
[u_lofi, z_lofi] = opt.Optimize(rand(n_z, 1));
time_lofioptimization = toc();
fprintf('Optimization finished in %.2f seconds\n', time_lofioptimization);

%% Visualize the ROM optimization results.

% Inspect the state solution.
u_lofi_reshape = reshape(u_lofi, n_yr, n_t);
Y_rom_1 = basis1.Decompress(u_lofi_reshape(1:r_1, :));
Y_rom_2 = basis2.Decompress(u_lofi_reshape(r_1 + 1:end, :));
Y_rom = [Y_rom_1; Y_rom_2];
solver.Animate_Solution(Y_rom);             % ROM state with ROM controller

% Inspect the control solution.
Q_rom = reshape(z_lofi, n_q, n_t - 1).^2;
figure;
semilogy(t(2:end), Q_rom);
title('Optimal controls (optimized with an OpInf ROM surrogate)');

%% Visualize the FOM with the ROM optimization results.

% Solve the high-fidelity model with the inferred controls.
disp('Final high-fidelity solve');
pp = spline(t(2:end), Q_rom);
controller = @(tt) ppval(pp, tt);
Y_hifi = solver.State_Solve(controller, t).NodalSolution;
solver.Animate_Solution(Y_hifi);            % FOM state with ROM controller

save('OptimizationSolution.mat', "solver", "Y_hifi", "Y_rom", "t", "Q_rom", "n_q");

%% Load and visualize results later.
% load('OptimizationSolution.mat', "solver", "Y_hifi", "Y_rom", "t", "Q_rom", "n_q");
% figure;
% plot(t(2:end), Q_rom);
% title('Optimal controls (optimized with an OpInf ROM surrogate)');
% solver.Animate_Solution(Y_rom);   % ROM state with ROM controller
% solver.Animate_Solution(Y_hifi);  % FOM state with ROM controller
