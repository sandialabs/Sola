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
data_interface = MD_Data_Interface_Transient_Test_Problem();
data_interface.Load_Data();
opt_prob_interface = MD_Opt_Prob_Interface_Sabl(opt, data_interface);
alpha_u = 1.e-2;
alpha_z = 1.e-12;
num_sing_vals = 90;
u_prior_interface = MD_Transient_Elliptic_u_Prior_Interface_Transient_Test_Problem(alpha_u, transient_prior_cov, opt, num_sing_vals);
z_prior_interface = MD_Elliptic_z_Prior_Interface_Transient_Test_Problem(alpha_z, opt);

%%
u_in = randn(n_y * n_t, 1);
u_out = u_prior_interface.Apply_W_u_Inverse_Factor(u_in);
tmp = reshape(u_in, n_y, n_t);
r = length(u_prior_interface.sing_vals);
test = sqrt(alpha_u) * u_prior_interface.sing_vecs_output * diag(u_prior_interface.sing_vals) * tmp(1:r, :) * diag(sqrt(transient_prior_cov.evals)) * transient_prior_cov.M_t_inv_evecs';

local_error = norm(u_out - test(:));
error = [error, local_error];

%%
F = u_prior_interface.Apply_W_u_Inverse_Factor(eye(n_y * n_t));
u_in = randn(n_y * n_t, 1);
u_out = u_prior_interface.Apply_W_u_Inverse(u_in);

local_error = norm(u_out - F * (F' * u_in));
error = [error, local_error];

%%
u_in = randn(n_y * n_t, 1);
scalar = rand;
u_out = u_prior_interface.Apply_W_u_Plus_scalar_M_u_Inverse_Factor(u_in, scalar);
aleph = kron(transient_prior_cov.evals, [u_prior_interface.sing_vals.^2; zeros(n_y - r, 1)]);
aleph = aleph ./ (1 + alpha_u * scalar * aleph);
test  = sqrt(alpha_u) * kron(transient_prior_cov.M_t_inv_evecs, [u_prior_interface.sing_vecs_output, zeros(n_y, n_y - r)]) * diag(sqrt(aleph)) * u_in;

local_error = norm(u_out - test);
error = [error, local_error];

%%
u_in = randn(n_y * n_t, 1);
scalar = rand;
F = u_prior_interface.Apply_W_u_Plus_scalar_M_u_Inverse_Factor(eye(n_y * n_t), scalar);
u_out = u_prior_interface.Apply_W_u_Plus_scalar_M_u_Inverse(u_in, scalar);

local_error = norm(u_out - F * (F' * u_in));
error = [error, local_error];

%%

if max(error) > 1.e-13
    disp('Error in model_discrepancy/Transient_Test_Problem/Driver_Unit_Test_2');
end
