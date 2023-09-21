clear
close all
clc

h = .1;
mesh = PDE_Meshing(h);
x = mesh.x;
y = mesh.y;
M = mesh.M;
diff_react_lofi = Diff_React_Lofi(mesh);

reg_coeff = 1.e-2;
diff_react_opt = Diff_React_Opt(diff_react_lofi,reg_coeff);

z0 = randn(diff_react_opt.m,1);
diff_react_opt.Finite_Difference_Gradient_Check(z0);
diff_react_opt.Finite_Difference_Hessian_Check(z0);