clear
close all
clc
addpath(genpath('../../src'))
rng(3524513)

load('Assembled_Operators.mat')

fd_check = false;
reg_coeff = 1.e-7;
adv_diff_opt = Adv_Diff_Opt(adv_diff,reg_coeff);

m = length(pde_meshing.x);
n = adv_diff_opt.n;
if fd_check
    z0 = 100*rand(n,1);
    adv_diff_opt.Finite_Difference_Gradient_Check(z0);
    adv_diff_opt.Finite_Difference_Hessian_Check(z0);
end

z0 = 10*ones(n,1);
[u_lofi,z_lofi] = adv_diff_opt.Optimize(z0);

T = adv_diff_opt.T;
u_hifi = nonlinear_adv_diff.State_Solve(adv_diff_opt.Map_z_vec_to_mesh(z_lofi));

name = 'Low-fidelity state';
pde_meshing.Plot_Field(u_lofi,name)

name = 'Low-fidelity control';
pde_meshing.Plot_Field(adv_diff_opt.Map_z_vec_to_mesh(z_lofi),name)
xlim(adv_diff_opt.control_xlim)
ylim(adv_diff_opt.control_ylim)

name = 'High-fidelity state';
pde_meshing.Plot_Field(u_hifi,name)

I = find(diag(adv_diff_opt.P_target)~=0);
cmin = min([u_lofi(I);T(I)]);
cmax = max([u_lofi(I);T(I)]);

name = 'Low-fidelity state';
pde_meshing.Plot_Field(u_lofi,name)
xlim(adv_diff_opt.target_xlim)
ylim(adv_diff_opt.target_ylim)
caxis([cmin,cmax])

name = 'Target state';
pde_meshing.Plot_Field(T,name)
xlim(adv_diff_opt.target_xlim)
ylim(adv_diff_opt.target_ylim)
caxis([cmin,cmax])

Z = z_lofi;
D = u_hifi-u_lofi;
save('Optimization_Results.mat','adv_diff','nonlinear_adv_diff','adv_diff_opt','z_lofi','u_lofi','Z','D')