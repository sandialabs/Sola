clear;
close all;
addpath(genpath('../../../src'));
rng(132253);

n_y = 100;
n_t = 51;
T = 1;
n_z = n_y;
obj = Adv_Diff_Objective(n_y, n_z, T, n_t);
con_hifi = Adv_Diff_Constraint(n_y, n_z, T, n_t);
con = Diff_Constraint(n_y, n_z, T, n_t);
opt = Reduced_Space_Optimization(obj, con);

data_interface = MD_Data_Interface_Transient_Test_Problem();
data_interface.Load_Data();

hyperparams = MD_Hyperparameters_Transient_Test_Problem(data_interface,n_y);
transient_prior_cov = MD_Transient_Prior_Covariance_Hyperparameters_Sabl(T, n_t, n_y, hyperparams);
u_prior_interface = MD_Transient_Laplacian_u_Prior_Interface(con_hifi.S,con_hifi.M,hyperparams,transient_prior_cov);
z_prior_interface = MD_Laplacian_z_Prior_Interface(con_hifi.S,con_hifi.M,hyperparams);


