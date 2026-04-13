%% Clear workspace and add path.
clear;
close all;
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

evalc('pcon.Parameterized_Finite_Difference_Constraint_Check(u, z, theta)');

%% Continuation algorithm
z0 = rand(2, 1);
evalc('[~, z_nom] = opt.Optimize(z0)');

sen_op = Euclidean_Sensitivity_Operators_Sabl(obj, pcon);
qn_prec = Quasi_Newton_Preconditioner();
sen = Pseudo_Time_Continuation(z_nom, sen_op, qn_prec);

evalc('sen_op.Finite_Difference_Gradient_Check(z, theta)');
evalc('sen_op.Finite_Difference_Hessian_Check(z, theta)');
evalc('sen_op.Finite_Difference_B_Check(z, theta)');

theta_pert = rand(4, 1);
opt.con.theta_current = theta_pert;
evalc('[~, z_pert] = opt.Optimize(z_nom)');

N = 200;
theta_traj = Euclidean_Auxillary_Parameter_Trajectory(N, theta_nom, theta_pert);
evalc('[z_k_fe, grad_k_fe] = sen.Pseudo_Time_Continuation_Forward_Euler(theta_traj)');
evalc('[z_k_me, grad_k_me] = sen.Pseudo_Time_Continuation_Modified_Euler(theta_traj)');

error = max(norm(z_pert - z_k_fe(:, end)), norm(z_pert - z_k_me(:, end)));

if error > 1.e-6
    fprintf(2, '\nparametric_sensitivities/Example_1 failed.\n');
else
    fprintf(1, '\nparametric_sensitivities/Example_1 passed.\n');
end
