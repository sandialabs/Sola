clear;
close all;
clc;
addpath(genpath('../../src'));

%% Define the optimization problem.
m = 200;
N = 51;
T = 1;
num_space_control_nodes = 10;
n = num_space_control_nodes * (N - 1);
obj = Adv_Diff_Gaussian_Source_Objective(m, n, T, N, num_space_control_nodes);
con = Adv_Diff_Gaussian_Source_Constraint(m, n, T, N, num_space_control_nodes);
opt = Reduced_Space_Optimization(obj, con);

%% Do finite difference checks with a random control.
z0 = randn(n, 1);
opt.Finite_Difference_Gradient_Check(z0);
opt.Finite_Difference_Hessian_Check(z0);
