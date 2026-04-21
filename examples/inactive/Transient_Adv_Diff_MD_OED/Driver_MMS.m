%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Method of Manufactured Solutions

clear;
close all;
clc;
addpath(genpath('../../src'));

plot_solution = true;

%% Set up the optimization problem.
n_y = 200;
n_t = 51;
T = 1;
n_z = n_t * n_y;
con = Adv_Diff_Constraint(n_y, n_z, T, n_t);

%% Define a custom control.
x = con.x;
t = con.t_mesh;
z0 = zeros(n_y, n_t);
for k = 1:n_t
    tk = t(k);
    z0(:, k) = cos(2 * pi * x) + 4 * pi^2 * tk * cos(2 * pi * x) - 2 * pi * tk * sin(2 * pi * x);
end
z0 = z0(:);

%% Solve the state equation with the custom control.
u = con.State_Solve(z0);
u_reshaped = reshape(u, n_y, n_t);

%% Compute the true state solution with the custom control.
u_true = zeros(n_y, n_t);
for k = 1:n_t
    tk = t(k);
    u_true(:, k) = tk * cos(2 * pi * x);
end

%% Plot the true and computed state solution.
if plot_solution
    figure;
    for k = 1:n_t
        plot(x, u_reshaped(:, k), '-', x, u_true(:, k), '--', 'LineWidth', 3);
        ylim([-1 1]);
        pause(.05);
    end
end

% Print the error of the computed solution.
disp(norm(u - u_true(:)) / norm(u_true(:)));
