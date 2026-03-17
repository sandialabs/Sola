%%
clear;
close all;
clc;
addpath(genpath('../../src'));
rng(123413);

m = 200;
x = linspace(0, 1, m)';
T = 20 * (x + .5) .* (1.3 - x);

z0 = rand(m, 1);
[u, z, D] = Optimal_Solution_Map(T, z0);

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

n_train = 50;
E = (1.e-2) * S + M;
tmp = linsolve(E, randn(m, n_train));
scaling = .4 * mean(abs(T)) / mean(abs(tmp(:)));
T_samples = T + scaling * tmp;

z_samples = zeros(m, n_train);
D_samples = zeros(m, m, n_train);
grad_norm = zeros(n_train, 1);
for k = 1:n_train
    [~, z_samples(:, k), D_samples(:, :, k), grad_norm(k)] = Optimal_Solution_Map(T_samples(:, k), z);
end
I = find(grad_norm < 1.e-7);
T_samples = T_samples(:, I) - T;
z_samples = z_samples(:, I) - z;
D_samples = D_samples(:, :, I);
n_train = size(T_samples, 2);

n_test = 20;
tmp = linsolve(E, randn(m, n_test));
scaling = .3 * mean(abs(T)) / mean(abs(tmp(:)));
T_test = T + scaling * tmp;
z_test = zeros(m, n_test);
grad_norm = zeros(n_test, 1);
for k = 1:n_test
    [~, z_test(:, k), ~, grad_norm(k)] = Optimal_Solution_Map(T_test(:, k), z);
end
I = find(grad_norm < 1.e-7);
T_test = T_test(:, I);
z_test = z_test(:, I);
n_test = size(T_test, 2);

%%
[U, Sigma, V] = svd(mean(D_samples, 3));
r = 32;
U_hdsa = U(:, 1:r);
V_hdsa = V(:, 1:r);
T_hdsa_samples = V_hdsa' * T_samples;
z_hdsa_samples = U_hdsa' * z_samples;

r1 = size(T_hdsa_samples, 1);
r2 = size(z_hdsa_samples, 1);
L_hdsa = zeros(r2, r1);
for k = 1:r2
    L_hdsa(k, :) = linsolve(T_hdsa_samples', z_hdsa_samples(k, :)');
end

%%
[U, Sigma, V] = svd(z_samples);
r = 32;
U_pod = U(:, 1:r);

[U, Sigma, V] = svd(T_samples);
r = 32;
V_pod = U(:, 1:r);

T_pod_samples = V_pod' * T_samples;
z_pod_samples = U_pod' * z_samples;

r1 = size(T_pod_samples, 1);
r2 = size(z_pod_samples, 1);
L_pod = zeros(r2, r1);
for k = 1:r2
    L_pod(k, :) = linsolve(T_pod_samples', z_pod_samples(k, :)');
end

%%
z_hdsa = zeros(m, n_test);
z_pod = zeros(m, n_test);
for k = 1:n_test
    z_hdsa(:, k) = z + U_hdsa * L_hdsa * V_hdsa' * (T_test(:, k) - T);
    z_pod(:, k) = z + U_pod * L_pod * V_pod' * (T_test(:, k) - T);
end

hdsa_error = zeros(n_test, 1);
pod_error = zeros(n_test, 1);
for k = 1:n_test
    hdsa_error(k) = sqrt((z_hdsa(:, k) - z_test(:, k))' * M * (z_hdsa(:, k) - z_test(:, k))) / sqrt(z_test(:, k)' * M * z_test(:, k));
    pod_error(k) = sqrt((z_pod(:, k) - z_test(:, k))' * M * (z_pod(:, k) - z_test(:, k))) / sqrt(z_test(:, k)' * M * z_test(:, k));
end

k = 1;
figure;
hold on;
plot(x, T_samples(:, 1), 'LineWidth', 3, 'color', [.9, .9, .9]);
plot(x, T_test(:, k), 'LineWidth', 3, 'color', 'magenta');
plot(x, T_samples, 'LineWidth', 3, 'color', [.9, .9, .9]);
plot(x, T_test(:, k), 'LineWidth', 3, 'color', 'magenta');
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$T$', 'Interpreter', 'latex');
legend({'Training Data', 'Testing Data'}, 'location', 'best', 'Interpreter', 'latex');
set(gca, 'FontSize', 24);
set(gcf, 'Color', 'White');

figure;
hold on;
plot(x, z_samples(:, 1), 'LineWidth', 3, 'color', [.9, .9, .9]);
plot(x, z_test(:, k), 'LineWidth', 3, 'color', 'magenta');
plot(x, z_samples, 'LineWidth', 3, 'color', [.9, .9, .9]);
plot(x, z_test(:, k), 'LineWidth', 3, 'color', 'magenta');
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$z$', 'Interpreter', 'latex');
legend({'Training Data', 'Testing Data'}, 'location', 'best', 'Interpreter', 'latex');
set(gca, 'FontSize', 24);
set(gcf, 'Color', 'White');

figure;
hold on;
plot(x, z_test(:, k), 'LineWidth', 3, 'color', 'magenta');
plot(x, z_pod(:, k), 'LineWidth', 3, 'color', 'cyan');
plot(x, z_hdsa(:, k), 'LineWidth', 3, 'color', 'red');
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$z$', 'Interpreter', 'latex');
legend({'Testing Data', 'POD Prediction', 'HDSA Prediction'}, 'location', 'best', 'Interpreter', 'latex');
set(gca, 'FontSize', 24);
set(gcf, 'Color', 'White');
