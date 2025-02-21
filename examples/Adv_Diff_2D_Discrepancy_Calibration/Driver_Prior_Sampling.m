%% Set up
clear;
close all;
clc;
addpath(genpath('../../src'));
load Optimization_Results.mat;
load Assembled_Operators.mat;
rng(2451423);

x = pde_meshing.x;
y = pde_meshing.y;
m = length(x);

con = Diff_Constraint(pde_meshing, diff_coeff);
obj = Diff_Objective(con, reg_coeff);
opt = Reduced_Space_Optimization(obj, con);

data_interface = MD_Data_Interface_Diff();
data_centering = true;
hyperparams = MD_Hyperparameters_Diff(data_interface, x, y, data_centering);

u_prior_interface = MD_Numeric_Laplacian_u_Prior_Interface(pde_meshing.S, pde_meshing.M, hyperparams);
z_prior_interface = MD_Numeric_Laplacian_z_Prior_Interface(pde_meshing.S, pde_meshing.M, hyperparams);

md_prior_sampling = MD_Prior_Sampling(data_interface, u_prior_interface, z_prior_interface);

num_prior_samples = 100;
num_perts = 10;
[delta_samples_z_opt, delta_samples_z_pert, z_pert] = md_prior_sampling.Prior_Discrepancy_Samples_for_Visualization(num_prior_samples, num_perts);

name = 'Discrepancy sample 1 at z_{opt}';
pde_meshing.Plot_Field(delta_samples_z_opt(:, 1), name);

name = 'Discrepancy sample 1 at pertubed z';
pde_meshing.Plot_Field(delta_samples_z_pert{1}(:, 1), name);

name = 'Perturbed z';
pde_meshing.Plot_Field(z_pert(:, 1), name);
