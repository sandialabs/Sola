%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
close all;
clc;
addpath(genpath('../../src'));

% Problem parameters
m = 200;             % number of spatial nodes
k0 = 1;              % reference permeability
alpha = 0.5;        % pressure-permeability sensitivity
p0 = 0;              % reference and boundary pressure
viscosity = 1;       % fluid viscosity
reg_coeff = 1.e-6;   % injection/production regularization

% Extra high-fidelity localized pressure-leakoff nonlinearity.  The low-
% fidelity model intentionally omits this stabilizing pressure sink, which
% creates a sizeable nonlinear model discrepancy for optimization and OED.
hifi_leakoff_coeff = 30;
hifi_leakoff_center = 0.55;
hifi_leakoff_width = 0.25;

obj = Subsurface_Objective(m, reg_coeff, p0);
con_lofi = Subsurface_LoFi_Constraint(m, k0, alpha, p0, viscosity);
con_hifi = Subsurface_HiFi_Constraint(con_lofi, hifi_leakoff_coeff, hifi_leakoff_center, hifi_leakoff_width);

opt_lofi = Reduced_Space_Optimization(obj, con_lofi);
opt_hifi = Reduced_Space_Optimization(obj, con_hifi);

% Solve the optimal control problem for the low- and high-fidelity models.
z0 = zeros(m, 1);
[u_lofi, z_lofi] = opt_lofi.Optimize(z0);
[u_hifi, z_hifi] = opt_hifi.Optimize(z_lofi);

% Report the high-fidelity objective gap induced by using the low-fidelity
% control.  These diagnostics are useful before launching the article/OED
% drivers, which use the same quantities at Step 0.
Jhat_lofi_on_hifi = opt_hifi.Jhat(z_lofi);
Jhat_hifi = opt_hifi.Jhat(z_hifi);
fprintf('High-fidelity objective at low-fidelity control: %.6e\n', Jhat_lofi_on_hifi);
fprintf('High-fidelity objective at high-fidelity control: %.6e\n', Jhat_hifi);
fprintf('Relative objective gap: %.2f%%\n', 100 * (Jhat_lofi_on_hifi - Jhat_hifi) / max(abs(Jhat_hifi), eps));

x = con_lofi.x;
T = obj.T;
u = con_hifi.State_Solve(z_lofi);
ymax = 1.1 * max([u_lofi; u_hifi; u; T]);
ymin = 1.1 * min([u_lofi; u_hifi; u; T; p0]);

figure;
hold on;
plot(x, T, 'LineWidth', 3);
plot(x, u, 'LineWidth', 3);
plot(x, u_lofi, 'LineWidth', 3);
xlabel('$x$', 'Interpreter', 'latex');
ylabel('Pressure', 'Interpreter', 'latex');
ylim([ymin, ymax]);
title('For low-fidelity injection control', 'Interpreter', 'latex');
legend({'Target', '$p$', '$\tilde{p}$'}, 'location', 'south', 'Interpreter', 'latex');
set(gca, 'FontSize', 24);
set(gcf, 'Color', 'White');

figure;
hold on;
plot(x, T, 'LineWidth', 3);
plot(x, u_hifi, 'LineWidth', 3);
xlabel('$x$', 'Interpreter', 'latex');
ylabel('Pressure', 'Interpreter', 'latex');
ylim([ymin, ymax]);
title('For high-fidelity injection control', 'Interpreter', 'latex');
legend({'Target', '$p$'}, 'location', 'south', 'Interpreter', 'latex');
set(gca, 'FontSize', 24);
set(gcf, 'Color', 'White');

figure;
hold on;
plot(x, z_hifi, 'LineWidth', 3, 'color', [0.8500 0.3250 0.0980]);
plot(x, z_lofi, 'LineWidth', 3, 'color', [0.9290 0.6940 0.1250]);
xlabel('$x$', 'Interpreter', 'latex');
ylabel('Injection/production rate', 'Interpreter', 'latex');
legend({'$q$', '$\tilde{q}$'}, 'location', 'north', 'Interpreter', 'latex');
set(gca, 'FontSize', 24);
set(gcf, 'Color', 'White');

% Generate a small set of control samples and evaluate the model discrepancy
% on them.  These samples are also used by the article/OED setup to
% automatically estimate prior hyperparameters, so include both the optimum
% and smooth nearby perturbations.
control_scale = max(abs(z_lofi));
if control_scale == 0
    control_scale = 1;
end
Z = zeros(m, 4);
Z(:, 1) = z_lofi;
Z(:, 2) = z_lofi + 0.25 * control_scale * sin(pi * x);
Z(:, 3) = z_lofi + 0.50 * control_scale * x .* (1 - x);
Z(:, 4) = z_lofi - 0.25 * control_scale * sin(2 * pi * x);

D = Evaluate_Discrepancy(con_hifi, con_lofi, Z);

save('Optimization_Results.mat', 'm', 'k0', 'alpha', 'p0', 'viscosity', ...
    'hifi_leakoff_coeff', 'hifi_leakoff_center', 'hifi_leakoff_width', ...
    'reg_coeff', 'z_lofi', 'z_hifi', 'u_lofi', 'u_hifi', 'Z', 'D');
