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

u_hyperparam_interface = MD_u_Hyperparameter_Interface_Diff(x, y, data_centering);
u_prior_interface = MD_Numeric_Laplacian_u_Prior_Interface(pde_meshing.S, pde_meshing.M, data_interface, u_hyperparam_interface);

num_state_solves = 100;
z_hyperparam_interface = MD_z_Hyperparameter_Interface_Diff(num_state_solves, x, y, con);
z_prior_interface = MD_Numeric_Laplacian_z_Prior_Interface(pde_meshing.S, pde_meshing.M, data_interface, z_hyperparam_interface, u_prior_interface);

% z_hyperparam_interface.beta_z = (5) * z_hyperparam_interface.beta_z;
% z_prior_interface.Set_beta_z(z_hyperparam_interface.beta_z);
% z_hyperparam_interface.alpha_z = (1/4) * z_hyperparam_interface.alpha_z;
% z_prior_interface.Set_alpha_z(z_hyperparam_interface.alpha_z);

md_prior_sampling = MD_Prior_Sampling(data_interface, u_prior_interface, z_prior_interface);
num_prior_samples = 100;
md_prior_sampling.Generate_Prior_Discrepancy_Sample_Data(num_prior_samples);

md_prior_vis = MD_Prior_Visualization(md_prior_sampling);
% md_prior_vis.Visualization_for_Prior_Discrepancy_at_z_opt(1);
md_prior_vis.Visualization_for_Prior_Discrepancy_at_z_pert(1);

% figure(3)
% xlim([0,1.3])
% set(figure(3), 'Position', [100, 100, 2200, 600]); % [left, bottom, width, height]
