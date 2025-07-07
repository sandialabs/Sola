clear;
close all;
clc;
addpath(genpath('../../src'));

m = load('Optimization_Results.mat', 'm').m;
vel_coeff = load('Optimization_Results.mat', 'vel_coeff').vel_coeff;
reg_coeff = load('Optimization_Results.mat', 'reg_coeff').reg_coeff;
diff_coeff = load('Optimization_Results.mat', 'diff_coeff').diff_coeff;
n_r = size(diff_coeff, 1);

obj = Adv_Diff_Objective(m, reg_coeff);
cons_hifi = cell(n_r, 1);
cons_lofi = cell(n_r, 1);
for k = 1:n_r
    cons_hifi{k} = Adv_Diff(m, vel_coeff, diff_coeff(k));
    cons_lofi{k} = Diff(cons_hifi{k});
end
opt_lofi = Reduced_Space_Optimization_Under_Uncertainty(obj, cons_lofi);

md_ouu_data_interface = MD_Data_Interface_Diff();
u_hyperparam_interface = load('Hyperparameter_Interfaces.mat', 'u_hyperparam_interface').u_hyperparam_interface;
z_hyperparam_interface = load('Hyperparameter_Interfaces.mat', 'z_hyperparam_interface').z_hyperparam_interface;

S = cons_lofi{1}.S;
M = cons_lofi{1}.M;

us_prior_interface = MD_Numeric_Laplacian_u_Prior_Interface(S, M, md_ouu_data_interface, u_hyperparam_interface);
u_prior_interface = MD_OUU_u_Prior_Interface(us_prior_interface, md_ouu_data_interface);
z_prior_interface = MD_Numeric_Laplacian_z_Prior_Interface(S, M, md_ouu_data_interface, z_hyperparam_interface, u_prior_interface);
