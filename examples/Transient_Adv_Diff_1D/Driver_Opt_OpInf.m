clear;
close all;
clc;
addpath(genpath('../../src'));

%% Experiment parameters.

datafile = 'OpInf_Training_Data.mat';
regenerate_data = false;
residual_energy_threshold = 1e-10;
regularization_candidates = logspace(-10, 3, 40);
ddt_strategy = '6thOrder';
plot_basis_functions = true;
plot_training_data = true;
plot_optimization_solution = true;
solve_hifi_problem = false;

%% Generate or load training data.

if ~exist(datafile, 'file') || regenerate_data
    disp('Generating training data');
    tic();

    % Set up the optimization problem.
    n_y = 200;
    n_t = 51;
    T = 1;
    n_q = 10;
    n_z = n_q * (n_t - 1);
    con = Adv_Diff_Gaussian_Source_Constraint(n_y, n_z, T, n_t, n_q);
    t = con.t_mesh;

    % Choose a target within the range of the PDE.
    q = zeros(n_q, n_t - 1);
    for i = 1:n_q
        q(i, :) = 10 * atan(10 * t(2:end)') .* cos(2 * pi * (t(2:end)' - (i - 1) / n_q));
    end
    z_truth = reshape(q, n_z, 1);
    target = reshape(con.State_Solve(z_truth), n_y, n_t);
    target_time_mesh = t;

    % Solve the state equation for several random controls.
    num_solves = 5;
    spln_nodes = 4;
    % Z_train = 10 * randn(n_z, num_solves);
    Z_train = zeros(n_z, num_solves);
    U_train = zeros(n_y * n_t, num_solves);
    for k = 1:num_solves
        vals = [zeros(n_q, 1), 5 * randn(n_q, spln_nodes - 1)];
        pp = spline(linspace(t(1), t(end), spln_nodes), vals);
        Z_train(:, k) = reshape(ppval(pp, t(2:end)), n_z, 1);
        U_train(:, k) = con.State_Solve(Z_train(:, k));
    end
    time_trainingdata = toc();

    % Save the experimental parameters and the training data.
    save(datafile, "n_y", "n_z", "T", "n_t", "n_q", "num_solves", ...
         "Z_train", "U_train", "time_trainingdata", "target", "target_time_mesh");
else
    load(datafile);
end
disp(['Using ', num2str(num_solves), ' training trajectories']);

%% Set up objectives and learn POD basis.

obj_hifi = Adv_Diff_Gaussian_Source_Objective(n_y, n_z, T, n_t, n_q);
t = obj_hifi.t_mesh;
x = obj_hifi.x;

states = reshape(U_train, n_y, n_t * num_solves);
basis = POD_Basis(states, false, obj_hifi.M);
basis.Set_Reduced_Dimension_From_Residual_Energy(residual_energy_threshold);
obj_lofi = Reduced_Dynamic_Objective(obj_hifi, basis.V);
disp(['Selected reduced dimension r = ', num2str(basis.r)]);

if plot_basis_functions
    figure;
    plot(x, basis.V(:, 1));
    hold on;
    for j = 2:basis.r
        plot(x, basis.V(:, j));
    end
    title('POD basis functions');
end

% Check projection error.
states_lofi = basis.Compress(states);
states_projected = basis.Decompress(states_lofi);
state_projection_error = norm(states_projected - states) / norm(states);
disp(['Projection error of training states: ', num2str(state_projection_error)]);
Yhats = reshape(states_lofi, basis.r, n_t, num_solves);
Qs = reshape(Z_train, n_q, [], num_solves);

if plot_training_data
    for k = 1:num_solves
        fig = figure();
        fig.Position(3:4) = [830, 300];
        subplot(1, 2, 1);
        plot(t(2:end), Qs(:, :, k));
        title(['training controls, trajectory', num2str(k)]);
        subplot(1, 2, 2);
        plot(t, Yhats(:, :, k));
        title(['compressed training states, trajectory ', num2str(k)]);
    end
end

target_projection_error = zeros(1, n_t - 1);
for k = 2:n_t
    target = obj_hifi.Evaluate_Target(obj_hifi.t_mesh(k), obj_hifi.x);
    target_projected = basis.Project(target);
    target_projection_error(k) = norm(target - target_projected) / norm(target);
end
disp(['Average projection error of target state: ', num2str(mean(target_projection_error))]);

%% Initialize and calibrate Operator Inference constraint.

operators = {Linear_Operator(), Input_Operator()};
rom = OpInf_ROM_Constraint(basis.r, n_q, T, n_t, zeros(basis.r, 1), operators);

tic();
rom.Select_Regularization(Yhats, Qs, regularization_candidates, ddt_strategy);
time_opinfcalibration = toc();

% Validate the Operator Inference constraint by solving
% the ROM for each of the training controls.
for k = 1:num_solves
    u = reshape(U_train(:, k), n_y, n_t);
    rom.y0 = Yhats(:, 1, k);
    u_rom_k = basis.Decompress(rom.State_Solve2(Qs(:, :, k)));
    state_error = norm(u - u_rom_k) / norm(u);
    disp(['ROM reconstruction error for training set ', num2str(k), ': ', num2str(100 * state_error), '%']);
end

%% Solve the optimization problem.

opt = Reduced_Space_Optimization(obj_lofi, rom);
z0 = randn(n_z, 1);

tic();
[u_reduced, z_lofi] = opt.Optimize(z0);
time_lofioptimization = toc();

u_rom = basis.Decompress(reshape(u_reduced, basis.r, n_t));

if ~solve_hifi_problem
    if plot_optimization_solution
        figure;
        for k = 1:n_t
            target = obj_hifi.Evaluate_Target(t(k), x);
            plot(x, u_rom(:, k), '-', x, target, '--', 'LineWidth', 3);
            title('Optimal states and target function');
            legend({'State (ROM)', 'Target'});
            xlim([0 1]);
            ylim([-.75 2]);
            pause(.1);
        end

        figure;
        hold on;
        z_lofi_reshape = reshape(z_lofi, n_q, n_t - 1);
        for i = 1:n_q
            plot(t(2:end), z_lofi_reshape(i, :), ...
                 '-', 'LineWidth', 1, 'Color', "#0072BD");
            plot(t(2:end), 10 * atan(10 * t(2:end)) .* cos(2 * pi * (t(2:end) - (i - 1) / n_q)), ...
                 '--', 'LineWidth', 1, 'Color', "#EDB120");
        end
        title('Optimal controls');
    end
    return
end

%% For comparison, solve the optimization problem in high fidelity.

con_hifi = Adv_Diff_Gaussian_Source_Constraint(n_y, n_z, T, n_t, n_q);
opt_hifi = Reduced_Space_Optimization(obj_hifi, con_hifi);

tic();
[u_true, z_hifi] = opt_hifi.Optimize(z0);
time_hifioptimization = toc();

u_true = reshape(u_true, n_y, n_t);
diff_state = norm(u_rom - u_true) / norm(u_true);
diff_control = norm(z_hifi - z_lofi) / norm(z_hifi);

%% Compare the optimal state to the target.

x = obj_hifi.x;
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
        title('Optimal states and target function');
        legend({'State (ROM)', 'State (FOM)', 'Target'});
        xlim([0 1]);
        ylim([-.75 2]);
        pause(.1);
    end
end
if plot_optimization_solution
    figure;
    semilogy(t, err1, '-', t, err2, ':', 'LineWidth', 3);
    legend({'Target error (ROM)', 'Target error (FOM)'});
    title('Error compared to target function');

    figure;
    hold on;
    z_hifi_reshape = reshape(z_hifi, n_q, n_t - 1);
    z_lofi_reshape = reshape(z_lofi, n_q, n_t - 1);
    for i = 1:n_q
        plot(t(2:end), z_lofi_reshape(i, :), ...
             '-', 'LineWidth', 1, 'Color', "#0072BD");
        plot(t(2:end), z_hifi_reshape(i, :), ...
             ':', 'LineWidth', 1, 'Color', "#D95319");
        plot(t(2:end), 10 * atan(10 * t(2:end)) .* cos(2 * pi * (t(2:end) - (i - 1) / n_q)), ...
             '--', 'LineWidth', 1, 'Color', "#EDB120");
    end
    title('Optimal controls');
end

%% Final report.

disp('Error Report');
disp(['  Relative error of low-fi and high-fi opt state:', 9, num2str(100 * diff_state), '%']);
disp(['  Relative error of low-fi and high-fi opt control:', 9, num2str(100 * diff_control), '%']);
disp(['  Relative error of low-fi opt state and target state:', 9, num2str(100 * mean(err1(2:end))), '%']);
disp(['  Relative error of high-fi opt state and target state:', 9, num2str(100 * mean(err2(2:end))), '%']);
disp('Timing Report');
disp(['  Generate ', num2str(num_solves), ' training trajectories:', 9, num2str(time_trainingdata), 's']);
disp(['  OpInf ROM calibration:', 9, 9, num2str(time_opinfcalibration), 's']);
disp(['  OpInf-constrained optimization:', 9, num2str(time_lofioptimization), 's']);
disp(['  FOM-constrained optimization:', 9, 9, num2str(time_hifioptimization), 's']);
