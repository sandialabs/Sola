%% Clear workspace and add path.
clear;
close all;
clc;
addpath(genpath('../../../src/optimization/'));
addpath(genpath('../../../src/parametric_sensitivities/'));
rng(1342);

%% Instantiate the optimization problem.
constants = [7; 1; 4; 8; 8];
obj = Example_1_Objective(constants);
theta_nom = [1; 1; 1; 1];
pcon = Example_1_Constraint(theta_nom);
opt = Reduced_Space_Optimization(obj, pcon);

%% Finite difference tests
u = rand(3, 1);
z = rand(2, 1);
theta = rand(4, 1);
[diffs_z, jacobian_z_transpose_check, diffs_u, jacobian_u_transpose_check, diff_theta, solve_res] = pcon.Parameterized_Finite_Difference_Constraint_Check(u, z, theta);

%% Continuation algorithm
z0 = rand(2, 1);
[~, z_nom] = opt.Optimize(z0);

sen_op = Sensitivity_Operators_Sabl(obj, pcon);
qn_prec = Quasi_Newton_Preconditioner();
sen = Pseudo_Time_Continuation(z_nom, theta_nom, sen_op, qn_prec);

sen_op.Finite_Difference_Gradient_Check(z, theta);
sen_op.Finite_Difference_Hessian_Check(z, theta);
sen_op.Finite_Difference_B_Check(z, theta);

theta_pert = rand(4, 1);
opt.con.theta_current = theta_pert;
[~, z_pert] = opt.Optimize(z_nom);

N = 200;
[z_k_fe, grad_k_fe] = sen.Pseudo_Time_Continuation_Forward_Euler(theta_pert, N);
[z_k_me, grad_k_me] = sen.Pseudo_Time_Continuation_Modified_Euler(theta_pert, N);
