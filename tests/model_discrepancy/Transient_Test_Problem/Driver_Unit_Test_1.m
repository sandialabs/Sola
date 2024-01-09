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

beta_tu = 11;
beta_iu = 1002;
beta_td = 4;

transient_prior_cov = MD_Transient_Prior_Covariance_Sabl(beta_tu, beta_iu, beta_td, T, n_t, n_y);

error = 0;

%%
E_tu_inv = linsolve(transient_prior_cov.E_tu, eye(n_t));
E_td_inv = linsolve(transient_prior_cov.E_td, eye(n_t));

local_error = norm(transient_prior_cov.evecs' * E_td_inv * transient_prior_cov.evecs - eye(n_t));
error = [error, local_error];

local_error = norm(E_td_inv * transient_prior_cov.evecs * diag(transient_prior_cov.evals) - E_tu_inv * transient_prior_cov.evecs) / norm(E_tu_inv * transient_prior_cov.evecs);
error = [error, local_error];

local_error = norm(E_td_inv * transient_prior_cov.evecs * diag(transient_prior_cov.evals) * transient_prior_cov.evecs' * E_td_inv - E_tu_inv) / norm(E_tu_inv);
error = [error, local_error];

%%
data_interface = MD_Data_Interface_Transient_Test_Problem();
data_interface.Load_Data();
opt_prob_interface = MD_Opt_Prob_Interface_Sabl(opt, data_interface);
alpha_u = 1.e-2;
alpha_z = 1.e-12;
num_sing_vals = 100;
u_prior_interface = MD_Elliptic_u_Prior_Interface_Transient_Test_Problem(alpha_u, transient_prior_cov, opt, num_sing_vals);
z_prior_interface = MD_Elliptic_z_Prior_Interface_Transient_Test_Problem(alpha_z, opt);

%%

E_tu = transient_prior_cov.E_tu;
E_td = transient_prior_cov.E_td;

W_su_inv = u_prior_interface.Apply_E_u_Inverse_Transpose(eye(n_y));
W_su_inv = u_prior_interface.Apply_M_u(W_su_inv);
W_su_inv = u_prior_interface.Apply_E_u_Inverse(W_su_inv);

W_su = linsolve(W_su_inv, eye(size(W_su_inv, 1)));

W_sd = u_prior_interface.Apply_E_d(eye(n_y));
W_sd = u_prior_interface.Apply_M_u_Inverse(W_sd);
W_sd = u_prior_interface.Apply_E_d_Transpose(W_sd);

W_u = (1 / alpha_u) * kron(E_tu, W_su);
W_d = kron(E_td, W_sd);

%%
local_error = norm(W_su_inv - u_prior_interface.sing_vecs_output * diag(u_prior_interface.sing_vals.^2) * u_prior_interface.sing_vecs_output');
error = [error, local_error];

%%
u_in = randn(n_y * n_t, 1);
u_out = u_prior_interface.Apply_W_d(u_in);
local_error = norm(u_out - W_d * u_in);
error = [error, local_error];

%%
u_in = randn(n_y * n_t, 1);
u_out = u_prior_interface.Apply_W_u_Inverse(u_in);
local_error = norm(u_out - linsolve(W_u, u_in));
error = [error, local_error];

%%
u_in = randn(n_y * n_t, 1);
u_out = u_prior_interface.Apply_W_u_Inverse_Factor(u_in);
tmp = reshape(u_in, n_y, n_t);
r = length(u_prior_interface.sing_vals);
test = sqrt(alpha_u) * u_prior_interface.sing_vecs_output * diag(u_prior_interface.sing_vals) * tmp(1:r, :) * diag(sqrt(transient_prior_cov.evals)) * transient_prior_cov.E_td_inv_evecs';

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
u_out = u_prior_interface.Apply_W_u_Plus_scalar_W_d_Inverse(u_in, scalar);
test  = linsolve(W_u + scalar * W_d, u_in);

local_error = norm(u_out - test);
error = [error, local_error];

%%
u_in = randn(n_y * n_t, 1);
scalar = rand;
u_out = u_prior_interface.Apply_W_u_Plus_scalar_W_d_Inverse_Factor(u_in, scalar);
aleph = kron(transient_prior_cov.evals, [u_prior_interface.sing_vals.^2; zeros(n_y - r, 1)]);
aleph = aleph ./ (1 + alpha_u * scalar * aleph);
test  = sqrt(alpha_u) * kron(transient_prior_cov.E_td_inv_evecs, [u_prior_interface.sing_vecs_output, zeros(n_y, n_y - r)]) * diag(sqrt(aleph)) * u_in;

local_error = norm(u_out - test);
error = [error, local_error];
%%

if max(error) > 1.e-7
    disp('Error in model_discrepancy/Transient_Test_Problem/Driver_Unit_Test_1');
end
