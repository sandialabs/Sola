clear;
close all;
clc;
addpath(genpath('../../src'));
rng(1423);

T = 8;
N = 100;
obj_hifi = Mass_Spring_Objective_HiFi(T, N);
obj_lofi = Mass_Spring_Objective_LoFi(obj_hifi);
con_hifi = Mass_Spring_Coupled(T, N);
con_lofi = Mass_Spring_LoFi(con_hifi);
opt_hifi = Reduced_Space_Optimization(obj_hifi, con_hifi);
opt_lofi = Reduced_Space_Optimization(obj_lofi, con_lofi);

z0 = rand(N - 1, 1);
fd_check = false;
if fd_check
    opt_hifi.Finite_Difference_Gradient_Check(z0);
    opt_hifi.Finite_Difference_Hessian_Check(z0);

    opt_lofi.Finite_Difference_Gradient_Check(z0);
    opt_lofi.Finite_Difference_Hessian_Check(z0);
end

[u_lofi, z_tilde] = opt_lofi.Optimize(z0);
u_hifi = con_hifi.State_Solve(z_tilde);

t = con_hifi.t_mesh;
u_lofi_t = reshape(u_lofi, [], N)';
u_hifi_t = reshape(u_hifi, [], N)';

[~, z_star] = opt_hifi.Optimize(z_tilde);

figure;
plot(t, u_lofi_t(:, 1), t, u_hifi_t(:, 1));

figure;
plot(t(2:end), z_tilde, t(2:end), z_star);

Z = z_tilde;
tmp = con_hifi.State_Solve(Z);
tmp = reshape(tmp, 14, N)';
tmp = tmp(:, 1:4)';
D = tmp(:) - con_lofi.State_Solve(Z);

save('Optimization_Results.mat', 'z_tilde', 'z_star', 'u_lofi', 'Z', 'D');
