%% Clear workspace and add path.
clear;
close all;
clc;
addpath('../../../src/optimization/');
rng(1342);

%% Instantiate the optimization problem.
constants = [7; 1; 4; 8; 8];
obj = Example_1_Objective(constants);
con = Example_1_Constraint();
opt = Reduced_Space_Optimization(obj, con);

%% Run finite difference checks.
opt.verbose = false;
z0 = rand(2, 1) + 1;
opt.Finite_Difference_Gradient_Check(z0);
opt.Finite_Difference_Hessian_Check(z0);

%% Do the optimization.
[u, z] = opt.Optimize(z0);

%% Compare solution to the truth and report discrepancies.
% The optimal solution is u = (a1  a2  a3), z = (a4  a5).
u_sol = constants(1:3);
z_sol = constants(4:5);

err = 0;
err = max(err, norm(u_sol - u));
err = max(err, norm(z_sol - z));
if err > 1e-12
    disp('Error in example 1');
    disp(err);
    disp('');
    disp('Calculated state');
    disp(u);
    disp('Calculated control');
    disp(z);
    disp('Calculated objective value');
    [val, gradu, gradz] = obj.J(u, z);
    disp(val);
    disp('');
    disp('True state');
    disp(u_sol);
    disp('True control');
    disp(z_sol);
    disp('Objective value at true solution');
    [val, gradu, gradz] = obj.J(u_sol, z_sol);
    disp(val);
end
