%% Clean slate.
clear;
close all;
clc;
run('../../src/Set_Paths');

%% Define the optimization objective and the high-fidelity constraint.

m = 200;
reg_coeff = 10;

% Initialize the optimization objective.
obj = Adv_Diff_Objective(m, reg_coeff);

% Set up an advection-diffusion system with some given system parameters.
diff_coeff = 1;
vel_coeff = 3;
robin_coeff = 2;
con_hifi = Adv_Diff_Constraint(m, diff_coeff, vel_coeff, robin_coeff);

% Solve a constrained optimization problem with the true system.
opt_hifi = Reduced_Space_Optimization(obj, con_hifi);
z0 = rand(m, 1);
[u_hifi, z_hifi] = opt_hifi.Optimize(z0);

%% Define the low-fidelity constraint and solve the optimization problem.

con_lofi = Diff_Constraint(con_hifi);
opt_lofi = Reduced_Space_Optimization(obj, con_lofi);
[u_lofi, z_lofi] = opt_lofi.Optimize(z0);

%% Plot results.

% Extract a few things for plotting.
x = con_hifi.x;
T = obj.T;
u = con_hifi.State_Solve(z_lofi);

% Compare the controls from each problem.
figure;
set(gcf, 'Position', [100, 100, 600, 150]);
subplot(1, 2, 1);
plot(x, z_hifi, 'LineStyle', '--', 'LineWidth', 1.5);
hold on;
plot(x, z_lofi, 'LineStyle', '-.', 'LineWidth', 1.5);
title('Controls', 'FontSize', 14, 'Interpreter', 'latex');
ylim([0, 1.1 * max([z_hifi; z_lofi]);]);
legend({'$z^*$', '$\tilde{z}$'}, 'Interpreter', 'latex');

subplot(1, 2, 2);
plot(x, u_hifi, 'LineStyle', '--', 'LineWidth', 1.5);
hold on;
plot(x, u_lofi, 'LineStyle', '-.', 'LineWidth', 1.5);
plot(x, u, 'LineStyle', ':', 'LineWidth', 1.5)
plot(x, T, 'Color', 'k', 'LineStyle', '-', 'LineWidth', 1.5);
title('States', 'Fontsize', 14, 'Interpreter', 'latex');
ylim([0, 1.1 * max([u_lofi; u_hifi; u; T]);]);
legend({'$s(z^{*})$', '$\tilde{s}(\tilde{z})$', '$s(\tilde{z})$', '$\phi$'}, ...
    'Location', 'southeast', 'Interpreter', 'latex');

%% AI CODE
% Define Matplotlib default colors
C0 = [0.1216, 0.4667, 0.7059]; % Blue
C1 = [1.0000, 0.4980, 0.0549]; % Orange
C2 = [0.1725, 0.6275, 0.1725]; % Green

% Create the figure
figure;
set(gcf, 'Position', [100, 100, 600, 150]); % Set figure size

% Use tiledlayout for better control over subplot spacing
tiledlayout(1, 2, 'Padding', 'compact', 'TileSpacing', 'compact');

% First subplot: Controls
nexttile;
plot(x, z_hifi, 'LineStyle', '--', 'LineWidth', 1.5, 'Color', C0); % High-fidelity control
hold on;
plot(x, z_lofi, 'LineStyle', '-.', 'LineWidth', 1.5, 'Color', C1); % Low-fidelity control
title('Controls', 'FontSize', 14, 'Interpreter', 'latex');
xlabel('$x$', 'FontSize', 12, 'Interpreter', 'latex'); % Add x-axis label with LaTeX
ylabel('$z(x)$', 'FontSize', 12, 'Interpreter', 'latex'); % Add y-axis label with LaTeX
ylim([0, 1.1 * max([z_hifi; z_lofi])]);
legend({'$z^*$', '$\tilde{z}$'}, 'Interpreter', 'latex', 'Box', 'off'); % Remove legend outline

% Second subplot: States
nexttile;
plot(x, T, 'Color', 'k', 'LineStyle', '-', 'LineWidth', 1.5); % Plot T first
hold on;
plot(x, u_hifi, 'LineStyle', '--', 'LineWidth', 1.5, 'Color', C0); % High-fidelity state
plot(x, u_lofi, 'LineStyle', '-.', 'LineWidth', 1.5, 'Color', C1); % Low-fidelity state

% Plot u with solid line and markers every 50 steps
markerIndices = 25:50:length(x); % Indices for markers every 50 steps
plot(x, u, 'LineStyle', '-', 'LineWidth', 1.5, 'Color', C2); % Solid line
plot(x(markerIndices), u(markerIndices), 'o', 'Color', C2, 'MarkerFaceColor', C2); % Markers

title('States', 'FontSize', 14, 'Interpreter', 'latex');
xlabel('$x$', 'FontSize', 12, 'Interpreter', 'latex'); % Add x-axis label with LaTeX
ylabel('$u(x)$', 'FontSize', 12, 'Interpreter', 'latex'); % Add y-axis label with LaTeX
ylim([0, 1.1 * max([u_lofi; u_hifi; u; T])]);
legend({'$\phi$', '$s(z^{*})$', '$\tilde{s}(\tilde{z})$', '$s(\tilde{z})$'}, ...
    'Location', 'southeast', 'Interpreter', 'latex', 'Box', 'off'); % Remove legend outline

% Export data for plotting in Python.
save('optimization_plot_data.mat', 'x', 'z_hifi', 'z_lofi', 'u_hifi', 'u_lofi', 'u', 'T');


%% Save data for later.
% Discrepancy data: surrogate control and constant control.
Z = zeros(m, 2);
Z(:, 1) = z_lofi;
Z(:, 2) = ones(m, 1);

% Discrepancy data s(z_k) - tilde{s}(z_k).
D = zeros(m, 2);
for k = 1:size(D, 2)
    D(:, k) = con_hifi.State_Solve(Z(:, k)) - con_lofi.State_Solve(Z(:, k));
end

save('Optimization_Results.mat', 'm', 'diff_coeff', 'vel_coeff', 'robin_coeff', 'reg_coeff', 'z_lofi', 'z_hifi', 'u_lofi', 'Z', 'D');
