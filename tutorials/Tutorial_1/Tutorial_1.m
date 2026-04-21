%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Clear workspace and add the SABL source path.
clear;
close all;
clc;
rng(1342);                              % Random seed for reproducing results (optional).

%% Instantiate the objective and constraints.
alphas = [7; 1; 4; 8; 8];
obj = Tutorial_1_Objective(alphas);
con = Tutorial_1_Constraint();

%% Run finite difference checks.
n_u = 3;
n_z = 2;
u0 = rand(n_u, 1) + 2;
z0 = rand(n_z, 1) + 1;

obj.Finite_Difference_Gradient_Check(u0, z0);
obj.Finite_Difference_Hessian_Check(u0, z0);
con.Finite_Difference_Constraint_Check(u0, z0);

%% Instantiate the optimization problem, solve it, and report the results.
opt = Reduced_Space_Optimization(obj, con);
[u, z] = opt.Optimize(z0);

disp("Experiment 1: alphas = [" + num2str(alphas') + "]");
disp(" ");
disp("State:");
disp(u);
disp("Control:");
disp(z);
disp("Objective:");
disp(obj.J(u, z));

% Try again with the Gauss-Newton approximation for the Hessian.
opt.Gauss_Newton_Hess = true;
[u, z] = opt.Optimize(z0);

disp("Experiment 2: alphas = [" + num2str(alphas') + "]");
disp(" ");
disp("State:");
disp(u);
disp("Control:");
disp(z);
disp("Objective:");
disp(obj.J(u, z));

%% Instantiate and solve the problem for different alphas.
alphas = [1; 2; 3; 4; 5];
obj = Tutorial_1_Objective(alphas);
opt = Reduced_Space_Optimization(obj, con);

[u, z] = opt.Optimize(z0);

% Report the results.
disp("Experiment 3: alphas = [" + num2str(alphas') + "]");
disp(" ");
disp("OPTIMIZATION");
disp("------------");
disp("State:");
disp(u);
disp("Control:");
disp(z);
disp("Objective:");
disp(obj.J(u, z));

%% Compare the results to a gridsearch of the control space.
minval = 1e10;
winner = [];
for z1 = linspace(2, 6, 500)
    for z2 = linspace(2, 6, 500)
        ztemp = [z1; z2];
        val = obj.J(con.State_Solve(ztemp), ztemp);
        if val < minval
            minval = val;
            winner = ztemp;
        end
    end
end

disp("GRIDSEARCH");
disp("----------");
disp("State:");
disp(con.State_Solve(winner));
disp("Control:");
disp(winner);
disp("Objective:");
disp(minval);

%% Use automatic differentiation classes to solve the optimization problem.
alphas = [7; 1; 4; 8; 8];
obj_AD = Tutorial_1_Objective_AD(alphas, n_u, n_z);
con_AD = Tutorial_1_Constraint_AD(n_u, n_z);
opt_AD = Reduced_Space_Optimization(obj_AD, con_AD);

evalc("obj_AD.AD_Initialization()");
evalc("con_AD.AD_Initialization()");

[u, z] = opt_AD.Optimize(z0);

% Report the results.
disp("Experiment 4: alphas = [" + num2str(alphas') + "]");
disp(" ");
disp("State:");
disp(u);
disp("Control:");
disp(z);
disp("Objective:");
disp(obj_AD.J(u, z));

% If new alphas are used, reinstantiate the objective.
alphas = [5; 4; 3; 2; 1];
obj_AD = Tutorial_1_Objective_AD(alphas, n_u, n_z);
opt_AD = Reduced_Space_Optimization(obj_AD, con_AD);
evalc("obj_AD.AD_Initialization()");

[u, z] = opt_AD.Optimize(z0);

disp("Experiment 5: alphas = [" + num2str(alphas') + "]");
disp(" ");
disp("State:");
disp(u);
disp("Control:");
disp(z);
disp("Objective:");
disp(obj_AD.J(u, z));
