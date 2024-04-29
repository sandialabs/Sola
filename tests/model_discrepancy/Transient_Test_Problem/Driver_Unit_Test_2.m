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

beta_t = 11;
beta_i = 1002;

transient_prior_cov = MD_Transient_Prior_Covariance_Sabl(beta_t, beta_i, T, n_t, n_y);

error = 0;

%%
data_interface = MD_Data_Interface_Transient_Test_Problem();
data_interface.Load_Data();
opt_prob_interface = MD_Opt_Prob_Interface_Sabl(opt, data_interface);
alpha_u = 1.e-2;
alpha_z = 1.e-12;
num_sing_vals = 90;
u_prior_interface = MD_Transient_Elliptic_u_Prior_Interface_Transient_Test_Problem(alpha_u, transient_prior_cov, opt, num_sing_vals);
z_prior_interface = MD_Elliptic_z_Prior_Interface_Transient_Test_Problem(alpha_z, opt);

%%
num_samples = 1.e4;

u_samps = u_prior_interface.Sample_with_Covariance_W_u_Inverse(num_samples);
test = cov(u_samps');
W_u_inv = u_prior_interface.Apply_W_u_Inverse(eye(n_y * n_t));
sampling_local_error = norm(W_u_inv - test, 'fro') / norm(W_u_inv, 'fro');
sampling_error = sampling_local_error;

%%
scalar = rand;
u_samps = u_prior_interface.Sample_with_Covariance_W_u_Plus_scalar_M_u_Inverse(num_samples, scalar);
test = cov(u_samps');
W_u_Plus_scalar_M_u_inv = u_prior_interface.Apply_W_u_Plus_scalar_M_u_Inverse(eye(n_y * n_t), scalar);
sampling_local_error = norm(W_u_Plus_scalar_M_u_inv - test, 'fro') / norm(W_u_Plus_scalar_M_u_inv, 'fro');
sampling_error = [sampling_error, sampling_local_error];

%%
if max(sampling_error) > .03
    disp('Error in model_discrepancy/Transient_Test_Problem/Driver_Unit_Test_1');
end
