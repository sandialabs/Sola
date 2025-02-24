clear;
close all;
clc;
addpath(genpath('../../src'));

n_y = 100;
n_t = 51;
T = 1;
n_z = n_y;
obj = Adv_Diff_Objective(n_y, n_z, T, n_t);
con_hifi = Adv_Diff_Constraint(n_y, n_z, T, n_t);
con = Diff_Constraint(n_y, n_z, T, n_t);
opt = Reduced_Space_Optimization(obj, con);

z0 = randn(n_y, 1);
[u_lofi, z_lofi] = opt.Optimize(z0);

Z = zeros(n_y, 2);
Z(:, 1) = z_lofi;
Z(:, 2) = con.x;

D = zeros(n_y * n_t, 2);
for k = 1:2
    D(:, k) = con_hifi.State_Solve(Z(:, k)) - con.State_Solve(Z(:, k));
end

save('Optimization_Results.mat', 'u_lofi', 'z_lofi', 'Z', 'D');
