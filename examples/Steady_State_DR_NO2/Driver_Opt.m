%%
clear;
close all;
clc;
addpath(genpath('../../src'));

m = 200;
x = linspace(0, 1, m)';
T = 20 * (x + .5) .* (1.3 - x);

z0 = rand(m, 1);
[u, z, D] = Optimal_Solution_Map(T, z0);

figure;
hold on;
plot(x, T, 'LineWidth', 3);
plot(x, u, 'LineWidth', 3);
xlabel('$x$', 'Interpreter', 'latex');
ylabel('Concentration', 'Interpreter', 'latex');
legend({'Target', '$u$'}, 'location', 'south', 'Interpreter', 'latex');
set(gca, 'FontSize', 24);
set(gcf, 'Color', 'White');

figure;
hold on;
plot(x, z, 'LineWidth', 3, 'color', [0.9290 0.6940 0.1250]);
xlabel('$x$', 'Interpreter', 'latex');
ylabel('Source', 'Interpreter', 'latex');
set(gca, 'FontSize', 24);
set(gcf, 'Color', 'White');

%%
execute_test = false;

if execute_test
    T_pert = T + 5;
    [u_pert, z_pert] = Optimal_Solution_Map(T_pert, z);
    z_pert_approx = z + D * (T_pert - T);

    figure;
    hold on;
    plot(x, z_pert, 'LineWidth', 3);
    plot(x, z_pert_approx, 'LineWidth', 3);
    xlabel('$x$', 'Interpreter', 'latex');
    ylabel('Source', 'Interpreter', 'latex');
    legend({'Optimal', 'Approximation'}, 'location', 'south', 'Interpreter', 'latex');
    set(gca, 'FontSize', 24);
    set(gcf, 'Color', 'White');
end

%%
h = x(2) - x(1);
M = diag(4 * ones(1, m)) + diag(ones(1, m - 1), 1) + diag(ones(1, m - 1), -1);
M(1, 1) = .5 * M(1, 1);
M(end, end) = .5 * M(end, end);
M = (1 / 6) * h * M;

S = diag(2 * ones(1, m)) + (-1) * diag(ones(1, m - 1), 1) + (-1) * diag(ones(1, m - 1), -1);
S(1, 1) = .5 * S(1, 1);
S(end, end) = .5 * S(end, end);
S = (1 / h) * S;

n = 10;
E = (1.e-2) * S + M;
tmp = linsolve(E, randn(m, n));
scaling = .2 * mean(abs(T)) / mean(abs(tmp(:)));
T_samples = T + scaling * tmp;

z_samples = 0 * T_samples;
for k = 1:n
    [~, z_samples(:, k)] = Optimal_Solution_Map(T_samples(:, k), z);
end

[U, Sigma, V] = svd(D);
r = 20;
Ur = U(:, 1:r);
Sigmar = Sigma(1:r, 1:r);
Vr = V(:, 1:r);

Tr_samples = Vr' * T_samples;
zr_samples = Ur' * z_samples;

L = zeros(r, r);
for k = 1:r
    L(k, :) = linsolve(Tr_samples', zr_samples(k, :)');
end

z_approx_1 = z + D * (T_samples - T);
z_approx_2 = z + Ur * Sigmar * Vr' * (T_samples - T);
z_approx_3 = z + Ur * L * Vr' * (T_samples - T);

for k = 1:n
    figure;
    hold on;
    plot(x, z_samples(:, k), 'LineWidth', 3);
    plot(x, z_approx_1(:, k), '--', 'LineWidth', 3);
    plot(x, z_approx_2(:, k), '*', 'LineWidth', 3);
    plot(x, z_approx_3(:, k), 'LineWidth', 3);
    xlabel('$x$', 'Interpreter', 'latex');
    ylabel('Source', 'Interpreter', 'latex');
    legend({'Optimal', 'Post-Opt Approx', 'LR Post-Opt Approx', 'NO2'}, 'location', 'south', 'Interpreter', 'latex');
    set(gca, 'FontSize', 24);
    set(gcf, 'Color', 'White');
end
