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

if false
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

n = 50;
E = (1.e-2) * S + M;
tmp = linsolve(E, randn(m, n));
scaling = .4 * mean(abs(T)) / mean(abs(tmp(:)));
T_samples = T + scaling * tmp;

z_samples = zeros(m,n);
D_samples = zeros(m,m,n);
grad_norm = zeros(n,1);
for k = 1:n
    [~, z_samples(:, k),D_samples(:,:,k),grad_norm(k)] = Optimal_Solution_Map(T_samples(:, k), z);
end
I = find(grad_norm<1.e-7);
T_samples = T_samples(:,I) - T;
z_samples = z_samples(:,I) - z;
D_samples = D_samples(:,:,I);

%%
[U, Sigma, V] = svd(mean(D_samples,3));
r = 32;
Ur = U(:, 1:r);
Sigmar = Sigma(1:r, 1:r);
Vr = V(:, 1:r);
Tr_samples = Vr' * T_samples;
zr_samples = Ur' * z_samples;

%%
L = zeros(r, r);
for k = 1:r
    L(k, :) = linsolve(Tr_samples', zr_samples(k, :)');
end

%%

% Define the number of hidden layers and neurons in each layer
p = 0; 
neurons_per_layer = r; 

% Create a feedforward neural network
net = feedforwardnet(repmat(neurons_per_layer, 1, p)); % p layers with r neurons each
net.biasConnect = false(1, net.numLayers); % Disable bias for all layers
net.trainFcn = 'trainlm'; % Levenberg-Marquardt backpropagation

% Divide data into training, validation, and test sets
net.divideParam.trainRatio = 80/100; % 70% for training
net.divideParam.valRatio = 20/100;   % 15% for validation
net.divideParam.testRatio = 0/100;  % 15% for testing

% % Set training parameters
% net.trainParam.epochs = 1000;          % Maximum number of epochs
% net.trainParam.goal = 1e-5;            % Performance goal
% net.trainParam.min_grad = 1e-16;        % Minimum gradient
% net.trainParam.max_fail = 100;            % Maximum validation failures
% net.trainParam.time = 60;              % Maximum training time in seconds

% Train the network
[net, tr] = train(net, Tr_samples, zr_samples);

%%
tmp = linsolve(E, randn(m, 1));
scaling = .3 * mean(abs(T)) / mean(abs(tmp(:)));
T_test = T + scaling * tmp;

z_approx = z + Ur * L * Vr' * (T_test - T);
[~, z_test, ~] = Optimal_Solution_Map(T_test, z);
z_no2 = z + Ur*net(Vr'*(T_test - T));

figure;
hold on;
plot(x, T_samples(:,1), 'LineWidth', 3, 'color', [.9, .9, .9]);
plot(x, T_test, 'LineWidth', 3, 'color', 'magenta');
plot(x, T_samples, 'LineWidth', 3, 'color', [.9, .9, .9]);
plot(x, T_test, 'LineWidth', 3, 'color', 'magenta');
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$T$', 'Interpreter', 'latex');
legend({'Training Data', 'Testing Data'}, 'location', 'best', 'Interpreter', 'latex');
set(gca, 'FontSize', 24);
set(gcf, 'Color', 'White');

figure;
hold on;
plot(x, z_samples(:,1), 'LineWidth', 3, 'color', [.9, .9, .9]);
plot(x, z_test, 'LineWidth', 3, 'color', 'magenta');
plot(x, z_no2, 'LineWidth', 3, 'color', 'cyan');
plot(x, z_approx, 'LineWidth', 3, 'color', 'red');
plot(x, z_samples, 'LineWidth', 3, 'color', [.9, .9, .9]);
plot(x, z_test, 'LineWidth', 3, 'color', 'magenta');
plot(x, z_no2, 'LineWidth', 3, 'color', 'cyan');
plot(x, z_approx, 'LineWidth', 3, 'color', 'red');
xlabel('$x$', 'Interpreter', 'latex');
ylabel('$z$', 'Interpreter', 'latex');
legend({'Training Data', 'Testing Data', 'NO2 Prediction','Linear Approx'}, 'location', 'best', 'Interpreter', 'latex');
set(gca, 'FontSize', 24);
set(gcf, 'Color', 'White');

