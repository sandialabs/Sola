clear;
close all;
clc;
addpath(genpath('../../src'));

m = 200;
diff_coeff = 1;
vel_coeff = 1 / 2;
robin_coeff = 2;
reg_coeff = 10;
obj = Adv_Diff_Objective(m, reg_coeff);
con_hifi = Adv_Diff_Constraint(m, diff_coeff, vel_coeff, robin_coeff);
con_lofi = Diff_Constraint(con_hifi);
opt_hifi = Reduced_Space_Optimization(obj, con_hifi);
opt_lofi = Reduced_Space_Optimization(obj, con_lofi);
x = con_hifi.x;

z0 = rand(m, 1);
[u_hifi, z_hifi] = opt_hifi.Optimize(z0);
[u_lofi, z_lofi] = opt_lofi.Optimize(z0);

x = con_hifi.x;
T = obj.T;
u = con_hifi.State_Solve(z_lofi);
ymax = 1.1 * max([u_lofi; u_hifi; u; T]);

figure;
hold on;
plot(x, T, 'LineWidth', 3);
plot(x, u, 'LineWidth', 3);
plot(x, u_lofi, 'LineWidth', 3);
xlabel('$x$', 'Interpreter', 'latex');
ylabel('Concentration', 'Interpreter', 'latex');
ylim([0, ymax]);
title('For low-fidelity control', 'Interpreter', 'latex');
legend({'Target', '$u$', '$\tilde{u}$'}, 'location', 'south', 'Interpreter', 'latex');
set(gca, 'FontSize', 24);
set(gcf, 'Color', 'White');

figure;
hold on;
plot(x, T, 'LineWidth', 3);
plot(x, u_hifi, 'LineWidth', 3);
xlabel('$x$', 'Interpreter', 'latex');
ylabel('Concentration', 'Interpreter', 'latex');
ylim([0, ymax]);
title('For high-fidelity control', 'Interpreter', 'latex');
legend({'Target', '$u$'}, 'location', 'south', 'Interpreter', 'latex');
set(gca, 'FontSize', 24);
set(gcf, 'Color', 'White');

figure;
hold on;
plot(x, z_hifi, 'LineWidth', 3, 'color', [0.8500 0.3250 0.0980]);
plot(x, z_lofi, 'LineWidth', 3, 'color', [0.9290 0.6940 0.1250]);
xlabel('$x$', 'Interpreter', 'latex');
ylabel('Source', 'Interpreter', 'latex');
legend({'$z$', '$\tilde{z}$'}, 'location', 'north', 'Interpreter', 'latex');
set(gca, 'FontSize', 24);
set(gcf, 'Color', 'White');

Z = zeros(m, 2);
Z(:, 1) = z_lofi;
Z(:, 2) = ones(m, 1);

D = zeros(m, 2);
for k = 1:size(D, 2)
    D(:, k) = con_hifi.State_Solve(Z(:, k)) - con_lofi.State_Solve(Z(:, k));
end

save('Optimization_Results.mat', 'm', 'diff_coeff', 'vel_coeff', 'robin_coeff', 'reg_coeff', 'z_lofi', 'z_hifi', 'u_lofi', 'Z', 'D');
