clear;
close all;
clc;
addpath(genpath('../../src'));

m = 100;
theta = ones(m, 1);
con = Adv_Diff_Constraint(theta);
x = con.x;

jacobian_check = true;
mms_check = true;
finite_diff_check = true;

if jacobian_check
    u = randn(m, 1);
    z = randn(m, 1);
    theta = randn(m, 1);
    [diffs_z, jacobian_z_transpose_check, diffs_u, jacobian_u_transpose_check, diffs_theta, solve_res] = con.Parameterized_Finite_Difference_Constraint_Check(u, z, theta);
end

if finite_diff_check
    Generate_Obs_Data(con);
    likelihood = Adv_Diff_Likelihood_Model(m);
    prior = Adv_Diff_Prior_Model(con);
    bayes_inv = Bayesian_Inversion(likelihood, prior, con);
    z0 = randn(m, 1);
    diffs = bayes_inv.opt.Finite_Difference_Gradient_Check(z0);
    diffs = bayes_inv.opt.Finite_Difference_Hessian_Check(z0);
end

if mms_check
    con.MMS_Check();
end
