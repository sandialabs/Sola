%%
clear;
close all;
addpath(genpath('../../../src'));
rng(121234);

con = Synthetic_Test_Constraint();
con_hifi = Synthetic_Test_Hifi_Constraint();
obj = Synthetic_Test_Objective(con);

opt = Reduced_Space_Optimization(obj, con);
opt_hifi = Reduced_Space_Optimization(obj, con_hifi);
z0 = randn(con.n, 1);
[u, z] = opt.Optimize(z0);
[~, z_hifi] = opt_hifi.Optimize(z0);

u_hifi = con_hifi.State_Solve(z);

figure;
hold on;
plot(u, 'LineWidth', 3);
plot(u_hifi, 'LineWidth', 3);

figure;
hold on;
plot(z, 'LineWidth', 3);
plot(z_hifi, 'LineWidth', 3);

Z = zeros(con.n, 3);
Z(:, 1) = z;
Z(:, 2:3) = randn(con.n, 2);

D = con_hifi.State_Solve(Z) - con.State_Solve(Z);

save('Opt_Data.mat', 'z', 'u', 'Z', 'D', 'z_hifi');
