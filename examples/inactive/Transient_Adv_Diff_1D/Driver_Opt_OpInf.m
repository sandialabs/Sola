%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
solve_hifi_problem = true;

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
    spln_nodes = 6;
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

% Unpack the states and controls by training trajectory.
states = cell(num_solves);
controls = cell(num_solves);
for k = 1:num_solves
    states{k} = reshape(U_train(:, k), n_y, n_t);
    controls{k} = reshape(Z_train(:, k), n_q, n_t - 1);
end

% Learn a POD basis from the collection of all state snapshots.
states_all = horzcat(states{:});
basis = POD_Basis(states_all, false, obj_hifi.M);
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

% Check the projection error for each trajectory.
states_lofi = cell(num_solves);
for k = 1:num_solves
    states_lofi{k} = basis.Compress(states{k});
    Yk_proj = basis.Decompress(states_lofi{k});
    proj_err_k = norm(Yk_proj - states{k}) / norm(states{k});
    disp(['Projection error of training states for trajectory ', ...
          num2str(k), ': ', ...
          num2str(proj_err_k)]);
end

if plot_training_data
    for k = 1:num_solves
        fig = figure();
        fig.Position(3:4) = [830, 300];
        subplot(1, 2, 1);
        plot(t(2:end), controls{k});
        title(['training controls, trajectory', num2str(k)]);
        subplot(1, 2, 2);
        plot(t, states_lofi{k});
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

y0 = zeros(basis.r, 1);
operators = {Linear_Operator(), Input_Operator()};
rom = OpInf_ROM_Constraint(basis.r, n_q, T, n_t, y0, operators);

tic();
rom.Select_Regularization(states_lofi, controls, regularization_candidates, ddt_strategy);
time_opinfcalibration = toc();

% Validate the Operator Inference constraint by solving
% the ROM for each of the training controls.
for k = 1:num_solves
    Yk_data = states{k};
    rom.y0 = states_lofi{k}(:, 1);
    Yk_rom = basis.Decompress(rom.State_Solve2(controls{k}));
    state_error = norm(Yk_data - Yk_rom) / norm(Yk_data);
    disp(['ROM reconstruction error for trajectory ', num2str(k), ': ', num2str(100 * state_error), '%']);
end
rom.y0 = y0;

%% Solve the optimization problem.

opt = Reduced_Space_Optimization(obj_lofi, rom);
z0 = randn(n_z, 1);

tic();
[u_lofi, z_lofi] = opt.Optimize(z0);
time_lofioptimization = toc();

Y_rom = basis.Decompress(reshape(u_lofi, basis.r, n_t));
u_rom = reshape(Y_rom, [], 1);

if ~solve_hifi_problem
    if plot_optimization_solution
        figure;
        for k = 1:n_t
            target = obj_hifi.Evaluate_Target(t(k), x);
            plot(x, Y_rom(:, k), '-', x, target, '--', 'LineWidth', 3);
            title('Optimal states and target function');
            legend({'State (ROM)', 'Target'});
            xlim([0 1]);
            ylim([-.75 2]);
            pause(.1);
        end

        figure;
        hold on;
        Q_rom = reshape(z_lofi, n_q, n_t - 1);
        for i = 1:n_q
            plot(t(2:end), Q_rom(i, :), ...
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

diff_state = norm(u_rom - u_true) / norm(u_true);
diff_control = norm(z_lofi - z_hifi) / norm(z_hifi);
Y_true = reshape(u_true, n_y, n_t);

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
    err1(k) = norm(target - Y_rom(:, k)) / denom;
    err2(k) = norm(target - Y_true(:, k)) / denom;
    if plot_optimization_solution
        plot(x, Y_rom(:, k), '-', x, Y_true(:, k), ':', x, target, '--', 'LineWidth', 3);
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
    Q_hifi = reshape(z_hifi, n_q, n_t - 1);
    Q_rom = reshape(z_lofi, n_q, n_t - 1);
    for i = 1:n_q
        plot(t(2:end), Q_rom(i, :), ...
             '-', 'LineWidth', 1, 'Color', "#0072BD");
        plot(t(2:end), Q_hifi(i, :), ...
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
