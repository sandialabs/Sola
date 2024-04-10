clear;
close all;
clc;
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
E_t_inv = linsolve(transient_prior_cov.E_t, eye(n_t));
M_t_inv = linsolve(transient_prior_cov.M_t, eye(n_t));

local_error = norm(transient_prior_cov.evecs' * M_t_inv * transient_prior_cov.evecs - eye(n_t));
error = [error, local_error];

local_error = norm(M_t_inv * transient_prior_cov.evecs * diag(transient_prior_cov.evals) - E_t_inv * transient_prior_cov.evecs) / norm(E_t_inv * transient_prior_cov.evecs);
error = [error, local_error];

local_error = norm(M_t_inv * transient_prior_cov.evecs * diag(transient_prior_cov.evals) * transient_prior_cov.evecs' * M_t_inv - E_t_inv) / norm(E_t_inv);
error = [error, local_error];

%%
data_interface = MD_Data_Interface_Transient_Test_Problem();
data_interface.Load_Data();
opt_prob_interface = MD_Opt_Prob_Interface_Sabl(opt, data_interface);
alpha_u = 1.e-2;
alpha_z = 1.e-12;
num_sing_vals = 100;
u_prior_interface = MD_Transient_Elliptic_u_Prior_Interface_Transient_Test_Problem(alpha_u, transient_prior_cov, opt, num_sing_vals);
z_prior_interface = MD_Elliptic_z_Prior_Interface_Transient_Test_Problem(alpha_z, opt);

%%

E_t = transient_prior_cov.E_t;
M_t = transient_prior_cov.M_t;

W_su_inv = u_prior_interface.Apply_E_u_Inverse_Transpose(eye(n_y));
W_su_inv = u_prior_interface.Apply_Spatial_M_u(W_su_inv);
W_su_inv = u_prior_interface.Apply_E_u_Inverse(W_su_inv);

W_su = linsolve(W_su_inv, eye(size(W_su_inv, 1)));

W_u = (1 / alpha_u) * kron(E_t, W_su);
M_u = kron(M_t, u_prior_interface.M);

%%
local_error = norm(W_su_inv - u_prior_interface.sing_vecs_output * diag(u_prior_interface.sing_vals.^2) * u_prior_interface.sing_vecs_output');
error = [error, local_error];

%%
u_in = randn(n_y * n_t, 1);
u_out = u_prior_interface.Apply_M_u(u_in);
local_error = norm(u_out - M_u * u_in);
error = [error, local_error];

%%
u_in = randn(n_y * n_t, 1);
u_out = u_prior_interface.Apply_W_u_Inverse(u_in);
local_error = norm(u_out - linsolve(W_u, u_in));
error = [error, local_error];

%%
u_in = randn(n_y * n_t, 1);
scalar = rand;
u_out = u_prior_interface.Apply_W_u_Plus_scalar_M_u_Inverse(u_in, scalar);
test  = linsolve(W_u + scalar * M_u, u_in);

local_error = norm(u_out - test);
error = [error, local_error];

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

if max(error) > 5.e-7
    disp('Error in model_discrepancy/Transient_Test_Problem/Driver_Unit_Test_1');
end

if max(sampling_error) > .03
    disp('Error in model_discrepancy/Transient_Test_Problem/Driver_Unit_Test_1');
end
