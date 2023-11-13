clear;
close all;
clc;
addpath(genpath('../../../src'));

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

alpha_u = 1.e-2;
alpha_z = 1.e-12;
hdsa = Diff_HDSA(opt, alpha_u, alpha_z, hdsa_trans_cov);

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
disp(['Maximum error = ', num2str(max(error))]);
