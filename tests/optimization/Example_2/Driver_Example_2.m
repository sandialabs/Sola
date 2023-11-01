clear;
close all;
addpath('../../../src/optimization/');
rng(142534);

m = 2;
n = 2;
T = .1;
N = 10^2;
obj = Example_2_Objective(m, n, T, N);
con = Example_2_Constraint(m, n, T, N);
opt = Reduced_Space_Optimization(obj, con);
opt.verbose = false;
z0 = randn(2, 1);
opt.Finite_Difference_Gradient_Check(z0);
opt.Finite_Difference_Hessian_Check(z0);
[u, z] = opt.Optimize(z0);

% The optimal solution should be
% u(1:2:end) \approx exp(t)
% u(2:2:end) \approx exp(t)
% z \approx [1 ; 1];

%%
u_sol = load('Solution_Example_2.mat', 'u').u;
z_sol = load('Solution_Example_2.mat', 'z').z;

error = 0;
error = max(error, norm(u_sol - u));
error = max(error, norm(z_sol - z));
if error ~= 0
    disp('Error in example 2');
end

% save('Solution_Example_2.mat','u','z','obj')
