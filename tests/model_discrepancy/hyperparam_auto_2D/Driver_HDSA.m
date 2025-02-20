%% Set up
clear;
close all;
clc;
addpath(genpath('../../src'));
load Optimization_Results.mat;
rng(2451423);

x = adv_diff.pde_meshing.x;
y = adv_diff.pde_meshing.y;
M = adv_diff.pde_meshing.M;
S = adv_diff.pde_meshing.S;
m = length(x);

obj = Adv_Diff_Objective(adv_diff, reg_coeff);
con = Adv_Diff_Constraint(adv_diff);
opt = Reduced_Space_Optimization(obj, con);

data_interface = MD_Data_Interface_Adv_Diff();
data_interface.Load_Data();

hyperparams = MD_Hyperparameters_hyperparam_auto_2D(data_interface, x, y);
u_prior_interface = MD_Numeric_Laplacian_u_Prior_Interface(S, M, hyperparams);

alpha_z = 1.e-4;
z_prior_interface = MD_Elliptic_z_Prior_Interface_Adv_Diff(alpha_z, opt);
