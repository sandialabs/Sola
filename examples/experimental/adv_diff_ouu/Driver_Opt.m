%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
close all;
clc;
addpath(genpath('../../src'));

m = 200;
vel_coeff = 1 / 2;
reg_coeff = 10;
n_r = 3;
diff_coeff = [1, .99, .8];
obj = Adv_Diff_Objective(m, reg_coeff);
cons_hifi = cell(n_r, 1);
cons_lofi = cell(n_r, 1);
for k = 1:n_r
    cons_hifi{k} = Adv_Diff(m, vel_coeff, diff_coeff(k));
    cons_lofi{k} = Diff(cons_hifi{k});
end

opt_hifi = Reduced_Space_Optimization_Under_Uncertainty(obj, cons_hifi);
opt_lofi = Reduced_Space_Optimization_Under_Uncertainty(obj, cons_lofi);

z0 = randn(m, 1);
[u_lofi, z_lofi] = opt_lofi.Optimize(z0);
[u_hifi, z_hifi] = opt_hifi.Optimize(z_lofi);

Z = zeros(m, 2);
Z(:, 1) = z_lofi;
Z(:, 2) = ones(m, 1);
D = zeros(m, n_r, 2);
for j = 1:2
    for k = 1:n_r
        D(:, k, j) = cons_hifi{k}.State_Solve(Z(:, j)) - cons_lofi{k}.State_Solve(Z(:, j));
    end
end

save('Optimization_Results.mat', 'Z', 'D', 'z_lofi', 'z_hifi', 'u_lofi', 'u_hifi', 'diff_coeff', 'm', 'vel_coeff', 'reg_coeff');
