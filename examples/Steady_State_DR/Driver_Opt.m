clear;
close all;
clc;
addpath(genpath('../../src'));

m = 200;
diff_coeff = 1;
react_coeff = -1;
reg_coeff = 1.e-4;
obj = Diff_React_Objective(m, reg_coeff);
con_lofi = Diff_React_Constraint(m, diff_coeff, react_coeff);
con_hifi = Diff_React_HiFi_Constraint(con_lofi);

opt_lofi = Reduced_Space_Optimization(obj, con_lofi);
opt_hifi = Reduced_Space_Optimization(obj, con_hifi);

z0 = rand(m, 1);
[u_lofi, z_lofi] = opt_lofi.Optimize(z0);
[u_hifi, z_hifi] = opt_hifi.Optimize(z_lofi);

x = con_lofi.x;
T = obj.T;
u = con_hifi.State_Solve(z_lofi);
ymax = 1.1 * max([u_lofi; u_hifi; u; T]);
% plot(x, u, 'LineWidth', 3);
figure;
hold on;
plot(x, T, "k--", 'LineWidth', 3);
plot(x, u_lofi, "r-", 'LineWidth', 3);
xlabel('$x$', 'Interpreter', 'latex');
ylabel('Concentration');
ylim([5, 18]);
title('Low-fidelity State and Target');
legend({'Target', '$\tilde{S}(\tilde{z})$'}, 'location', 'southeast', 'Interpreter', 'latex');
set(gca, 'FontSize', 20);
saveas(gcf, "Lo-Fi State.svg");


figure;
hold on;
plot(x, T, "k--", 'LineWidth', 3);
plot(x, u, "r-", 'LineWidth', 3);
xlabel('$x$', 'Interpreter', 'latex');
ylabel('Concentration');
ylim([5, 18]);
title('High-fidelity State and Target');
legend({'Target', '$S(\tilde{z})$'}, 'location', 'southeast', 'Interpreter', 'latex');
set(gca, 'FontSize', 20);
saveas(gcf, "Hi-Fi State.svg");
% set(gcf, 'Color', 'White');

figure;
hold on;
plot(x, T, 'LineWidth', 3);
plot(x, u_hifi, 'LineWidth', 3);
xlabel('$x$', 'Interpreter', 'latex');
ylabel('Concentration', 'Interpreter', 'latex');
ylim([0, ymax]);
title('For high-fidelity control', 'Interpreter', 'latex');
set(gca, 'FontSize', 24);
legend({'Target', '$u$'}, 'location', 'south', 'Interpreter', 'latex');
set(gcf, 'Color', 'White');

figure;
hold on;
% plot(x, z_hifi, "k--", 'LineWidth', 3);
plot(x, z_lofi, "r-", 'LineWidth', 3);
xlabel('$x$', 'Interpreter', 'latex');
ylabel('Source');
legend({'Lo-Fi Sol.'}, 'location', 'northeast', 'Interpreter', 'latex');
title("Optimal Control Solutions")
set(gca, 'FontSize', 20);
saveas(gcf, "OC-Sols.svg");
set(gcf, 'Color', 'White');

Z = zeros(m, 2);
Z(:, 1) = z_lofi;
Z(:, 2) = 4.5 * max(z_lofi) * x .* (1 - x);

D = zeros(m, 2);
for k = 1:size(D, 2)
    D(:, k) = con_hifi.State_Solve(Z(:, k)) - con_lofi.State_Solve(Z(:, k));
end

save('Optimization_Results.mat', 'm', 'diff_coeff', 'react_coeff', 'reg_coeff', 'z_lofi', 'z_hifi', 'u_lofi', 'Z', 'D');
