%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
close all;
clc;
addpath(genpath('../../src'));

m = load('Optimization_Results.mat', 'm').m;
vel_coeff = load('Optimization_Results.mat', 'vel_coeff').vel_coeff;
reg_coeff = load('Optimization_Results.mat', 'reg_coeff').reg_coeff;
diff_coeff = load('Optimization_Results.mat', 'diff_coeff').diff_coeff;
n_r = size(diff_coeff, 2);

cons_hifi = cell(n_r, 1);
cons_lofi = cell(n_r, 1);
for k = 1:n_r
    cons_hifi{k} = Adv_Diff(m, vel_coeff, diff_coeff(k));
    cons_lofi{k} = Diff(cons_hifi{k});
end
x = cons_lofi{1}.x;
S = cons_lofi{1}.S;
M = cons_lofi{1}.M;

md_ouu_data_interface = MD_Data_Interface_Diff();
md_ouu_hyperparam_data_interface = MD_OUU_Hyperparam_Data_Interface(md_ouu_data_interface);

data_centering = false;

u_hyperparam_interface = MD_u_Hyperparameter_Interface_Diff(x, data_centering);
u_hyperparam_interface.Set_beta_u(0.08);
u_prior_interface = MD_Numeric_Laplacian_u_Prior_Interface(S, M, md_ouu_hyperparam_data_interface, u_hyperparam_interface);

num_state_solves = 100;
z_hyperparam_interface = MD_z_Hyperparameter_Interface_Diff(num_state_solves, x, cons_lofi);
z_prior_interface = MD_Numeric_Laplacian_z_Prior_Interface(S, M, md_ouu_hyperparam_data_interface, z_hyperparam_interface, u_prior_interface);
z_hyperparam_interface.Set_alpha_z(5.0);

save('Hyperparameter_Interfaces.mat', 'u_hyperparam_interface', 'z_hyperparam_interface');
