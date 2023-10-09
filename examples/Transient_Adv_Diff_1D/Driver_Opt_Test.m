clear
close all
clc
addpath(genpath('../../src'))

%% Define the optimization problem.
m = 200;
N = 51;
T = 1;
num_space_control_nodes = 10;
n = num_space_control_nodes*(N-1);
adv_diff = Adv_Diff_Gaussian_Source(m,n,T,N,num_space_control_nodes);

%% Do finite difference checks with a random control.
z0 = randn(n,1);
adv_diff.Finite_Difference_Gradient_Check(z0);
adv_diff.Finite_Difference_Hessian_Check(z0);
