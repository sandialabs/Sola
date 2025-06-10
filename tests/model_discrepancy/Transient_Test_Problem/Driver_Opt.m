clear;
close all;
clc;
addpath(genpath('../../../src'));

n_y = 50;
n_t = 11;
T = 1;
n_z = n_y;
obj = Adv_Diff_Objective(n_y, n_z, T, n_t);
obj.w(1) = 0;
obj.w(end) = obj.w(end - 1);
con_hifi = Adv_Diff_Constraint(n_y, n_z, T, n_t);
con = Diff_Constraint(n_y, n_z, T, n_t);
opt = Reduced_Space_Optimization(obj, con);

z0 = ones(n_y, 1);
[u_lofi, z_lofi] = opt.Optimize(z0);

Z = z_lofi;
D = con_hifi.State_Solve(Z) - con.State_Solve(Z);

save('Optimization_Results.mat', 'u_lofi', 'z_lofi', 'Z', 'D');
