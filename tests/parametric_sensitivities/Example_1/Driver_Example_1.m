%% Clear workspace and add path.
clear;
close all;
clc;
addpath('../../../src/optimization/');
addpath('../../../src/parametric_sensitivities/');
rng(1342);

%% Instantiate the optimization problem.
constants = [7; 1; 4; 8; 8];
obj = Example_1_Objective(constants);
theta_nom = [1; 1; 1; 1];
pcon = Example_1_Constraint(theta_nom);
opt = Reduced_Space_Optimization(obj, pcon);
psen_op = Parameteric_Sensitivity_Operators(obj, pcon);

%% Finite difference tests
u = rand(3, 1);
z = rand(2, 1);
theta = rand(4, 1);
[diffs_z, jacobian_z_transpose_check, diffs_u, jacobian_u_transpose_check, diff_theta, solve_res] = pcon.Parameterized_Finite_Difference_Constraint_Check(u, z, theta);
diffs_H = psen_op.Finite_Difference_Hessian_Check(z, theta);
diffs_B = psen_op.Finite_Difference_B_Check(z, theta);

%% Continuation algorithm
z0 = rand(2, 1);
[~, z_nom] = opt.Optimize(z0);

theta_pert = rand(4, 1);
opt.con.theta_current = theta_pert;
[~, z_pert] = opt.Optimize(z_nom);

N = 100;
% pseudo_time_con = Pseudo_Time_Continuation(obj,pcon,z_nom,theta_nom);
pseudo_time_con = Pseudo_Time_Continuation_BFGS_Example_1(obj, pcon, z_nom, theta_nom);
z_pert_con = pseudo_time_con.Forward_Euler_Continuation(theta_pert, N);
