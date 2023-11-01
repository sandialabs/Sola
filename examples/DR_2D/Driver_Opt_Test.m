clear;
close all;
clc;
addpath(genpath('../../src'));

h = .1;
mesh = PDE_Meshing(h);
x = mesh.x;
y = mesh.y;
M = mesh.M;
diff_react_lofi = Diff_React_Lofi(mesh);

reg_coeff = 1.e-2;
obj = Diff_React_Objective(diff_react_lofi, reg_coeff);
con = Diff_React_Constraint(diff_react_lofi);
opt = Reduced_Space_Optimization(obj, con);

z0 = randn(length(x), 1);
opt.Finite_Difference_Gradient_Check(z0);
opt.Finite_Difference_Hessian_Check(z0);
