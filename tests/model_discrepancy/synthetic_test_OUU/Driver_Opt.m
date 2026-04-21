%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Clear workspace and add path.
clear;
close all;
clc;
rng(1342);

%% Instantiate the optimization problem.
N = 30;
Xi = ones(3, 1) * linspace(0, 1, N);

m = 51;
obj = Synthetic_Test_OUU_Objective(m);
cons = cell(N, 1);
for k = 1:N
    cons{k} = Synthetic_Test_OUU_Constraint(Xi(:, k));
end
opt = Reduced_Space_Optimization_Under_Uncertainty(obj, cons);

%% Run finite difference checks.
% z0 = rand(m, 1) + 1;
% opt.Finite_Difference_Gradient_Check(z0);
% opt.Finite_Difference_Hessian_Check(z0);

%% Do the optimization.
z0 = rand(m, 1) + 1;
[u_opt, z_opt] = opt.Optimize(z0);

%% Generation hifi-data
cons_hifi = cell(N, 1);
for k = 1:N
    cons_hifi{k} = Synthetic_Test_OUU_HiFi_Constraint(Xi(:, k));
end
opt_hifi = Reduced_Space_Optimization_Under_Uncertainty(obj, cons_hifi);
% z0 = rand(m, 1) + 1;
% opt_hifi.Finite_Difference_Gradient_Check(z0);
% opt_hifi.Finite_Difference_Hessian_Check(z0);
[u_hifi, z_hifi] = opt_hifi.Optimize(z_opt);

x = obj.x;
Z = zeros(m, 2);
Z(:, 1) = z_opt;
Z(:, 2) = x + x.^2;

D = zeros(m, N, 2);
for i = 1:2
    for k = 1:N
        D(:, k, i) = cons_hifi{k}.State_Solve(Z(:, i)) - cons{k}.State_Solve(Z(:, i));
    end
end

save('Optimization_Results.mat', 'u_opt', 'z_opt', 'Z', 'D', 'Xi', 'z_hifi');
