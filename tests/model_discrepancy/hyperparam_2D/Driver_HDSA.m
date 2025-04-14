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

u_hyperparams = MD_u_Hyperparameter_Interface_hyperparam_2D(x, y);
u_prior_interface = MD_Numeric_Laplacian_u_Prior_Interface(S, M, data_interface, u_hyperparams);
