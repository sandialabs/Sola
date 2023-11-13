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

hdsa_trans_cov = HDSA_Sabl_Transient_Prior_Covariance(beta_tu, beta_iu, beta_td, T, n_t);

error = 0;

%%
E_tu_inv = linsolve(hdsa_trans_cov.E_tu, eye(n_t));
E_td_inv = linsolve(hdsa_trans_cov.E_td, eye(n_t));

local_error = norm(hdsa_trans_cov.evecs' * E_td_inv * hdsa_trans_cov.evecs - eye(n_t));
error = [error, local_error];

local_error = norm(E_td_inv * hdsa_trans_cov.evecs * diag(hdsa_trans_cov.evals) - E_tu_inv * hdsa_trans_cov.evecs) / norm(E_tu_inv * hdsa_trans_cov.evecs);
error = [error, local_error];

local_error = norm(E_td_inv * hdsa_trans_cov.evecs * diag(hdsa_trans_cov.evals) * hdsa_trans_cov.evecs' * E_td_inv - E_tu_inv) / norm(E_tu_inv);
error = [error, local_error];

%%
alpha_u = 1.e-2;
alpha_z = 1.e-12;
md_interface = Diff_HDSA(opt, alpha_u, alpha_z, hdsa_trans_cov);

%%

E_tu = hdsa_trans_cov.E_tu;
E_td = hdsa_trans_cov.E_td;

W_su_inv = md_interface.Apply_E_u_Inverse_Transpose(eye(n_y));
W_su_inv = md_interface.Apply_M_u(W_su_inv);
W_su_inv = md_interface.Apply_E_u_Inverse(W_su_inv);

W_su = linsolve(W_su_inv, eye(size(W_su_inv, 1)));

W_sd = md_interface.Apply_E_d(eye(n_y));
W_sd = md_interface.Apply_M_u_Inverse(W_sd);
W_sd = md_interface.Apply_E_d_Transpose(W_sd);

W_u = (1 / alpha_u) * kron(E_tu, W_su);
W_d = kron(E_td, W_sd);

%%
local_error = norm(W_su_inv - md_interface.sing_vecs_output * diag(md_interface.sing_vals.^2) * md_interface.sing_vecs_output');
error = [error, local_error];

%%
u_in = randn(n_y * n_t, 1);
u_out = md_interface.Apply_W_d(u_in);
local_error = norm(u_out - W_d * u_in);
error = [error, local_error];

%%
u_in = randn(n_y * n_t, 1);
u_out = md_interface.Apply_W_u_Inverse(u_in);
local_error = norm(u_out - linsolve(W_u, u_in));
error = [error, local_error];

%%
u_in = randn(n_y * n_t, 1);
u_out = md_interface.Apply_W_u_Inverse_Factor(u_in);
test = sqrt(alpha_u) * kron(hdsa_trans_cov.E_td_inv_evecs * diag(sqrt(hdsa_trans_cov.evals)), md_interface.sing_vecs_output * diag(md_interface.sing_vals)) * u_in;

local_error = norm(u_out - test);
error = [error, local_error];

%%
u_in = randn(n_y * n_t, 1);
scalar = rand;
u_out = md_interface.Apply_W_u_Plus_scalar_W_d_Inverse(u_in, scalar);
test  = linsolve(W_u + scalar * W_d, u_in);

local_error = norm(u_out - test);
error = [error, local_error];

%%
u_in = randn(n_y * n_t, 1);
scalar = rand;
u_out = md_interface.Apply_W_u_Plus_scalar_W_d_Inverse_Factor(u_in, scalar);
aleph = kron(hdsa_trans_cov.evals, md_interface.sing_vals.^2);
aleph = aleph ./ (1 + alpha_u * scalar * aleph);
test  = sqrt(alpha_u) * kron(hdsa_trans_cov.E_td_inv_evecs, md_interface.sing_vecs_output) * diag(sqrt(aleph)) * u_in;

local_error = norm(u_out - test);
error = [error, local_error];
%%

if max(error) > 1.e-7
    disp('Error in model_discrepancy/Transient_Test_Problem/Driver_Unit_Test');
end
