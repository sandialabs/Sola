clear;
close all;
clc;
addpath(genpath('../../src'));
rng(3524513);

load('Assembled_Operators.mat');

fd_check = false;
reg_coeff = 1.e-7;
obj = Adv_Diff_Objective(adv_diff, reg_coeff);
con = Adv_Diff_Constraint(adv_diff);
opt = Reduced_Space_Optimization(obj, con);

m = length(pde_meshing.x);
n = obj.n;
if fd_check
    z0 = 100 * rand(n, 1);
    opt.Finite_Difference_Gradient_Check(z0);
    opt.Finite_Difference_Hessian_Check(z0);
end

z0 = 10 * ones(n, 1);
[u_lofi, z_lofi] = opt.Optimize(z0);

T = obj.T;
u_hifi = nonlinear_adv_diff.State_Solve(con.Map_z_vec_to_mesh(z_lofi));

name = 'Low-fidelity state';
pde_meshing.Plot_Field(u_lofi, name);

name = 'Low-fidelity control';
pde_meshing.Plot_Field(con.Map_z_vec_to_mesh(z_lofi), name);
xlim(obj.control_xlim);
ylim(obj.control_ylim);

name = 'High-fidelity state';
pde_meshing.Plot_Field(u_hifi, name);

I = find(diag(obj.P_target) ~= 0);
cmin = min([u_lofi(I); T(I)]);
cmax = max([u_lofi(I); T(I)]);

name = 'Low-fidelity state';
pde_meshing.Plot_Field(u_lofi, name);
xlim(obj.target_xlim);
ylim(obj.target_ylim);
caxis([cmin, cmax]);

name = 'Target state';
pde_meshing.Plot_Field(T, name);
xlim(obj.target_xlim);
ylim(obj.target_ylim);
caxis([cmin, cmax]);

Z = z_lofi;
D = u_hifi - u_lofi;
save('Optimization_Results.mat', 'adv_diff', 'nonlinear_adv_diff', 'reg_coeff', 'z_lofi', 'u_lofi', 'Z', 'D');
