%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
close all;
clc;
addpath(genpath('../../src'));
rng(3524513);

load('Assembled_Operators.mat');

fd_check = false;
reg_coeff = 1.e-5;
con = Diff_Constraint(pde_meshing, diff_coeff);
obj = Diff_Objective(con, reg_coeff);
opt = Reduced_Space_Optimization(obj, con);
opt.max_cg_iter = 500;

m = length(pde_meshing.x);
if fd_check
    z0 = 100 * rand(m, 1);
    opt.Finite_Difference_Gradient_Check(z0);
    opt.Finite_Difference_Hessian_Check(z0);
end

z0 = 10 * ones(m, 1);
[u_lofi, z_lofi] = opt.Optimize(z0);

T = obj.T;
u_hifi = adv_diff.State_Solve(z_lofi);

name = 'Low-fidelity state';
pde_meshing.Plot_Field(u_lofi, name);
caxis([0, 15]);
set(gca, 'fontsize', 24);

name = 'Low-fidelity control';
pde_meshing.Plot_Field(z_lofi, name);
set(gca, 'fontsize', 24);

name = 'High-fidelity state';
pde_meshing.Plot_Field(u_hifi, name);
caxis([0, 15]);
set(gca, 'fontsize', 24);

name = 'Target state';
pde_meshing.Plot_Field(T, name);
caxis([0, 15]);
set(gca, 'fontsize', 24);

name = 'Model Discrepancy';
pde_meshing.Plot_Field(u_hifi - u_lofi, name);
set(gca, 'fontsize', 24);

z1 = z_lofi;
Z = z_lofi;
D = u_hifi - u_lofi;

save('Optimization_Results.mat', 'adv_diff', 'con', 'reg_coeff', 'z_lofi', 'u_lofi', 'Z', 'D');
