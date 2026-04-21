%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
close all;
clc;
addpath(genpath('../../src'));

%% Define the optimization problem.
n_y = 200;
n_t = 51;
T = 1;
num_space_control_nodes = 10;
n_z = num_space_control_nodes * (n_t - 1);
obj = Adv_Diff_Gaussian_Source_Objective(n_y, n_z, T, n_t, num_space_control_nodes);
con = Adv_Diff_Gaussian_Source_Constraint(n_y, n_z, T, n_t, num_space_control_nodes);
opt = Reduced_Space_Optimization(obj, con);

%% Do finite difference checks with a random control.
z0 = randn(n_z, 1);
opt.Finite_Difference_Gradient_Check(z0);
opt.Finite_Difference_Hessian_Check(z0);
