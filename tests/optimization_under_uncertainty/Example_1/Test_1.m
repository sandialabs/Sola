%% Clear workspace and add path.
clear;
close all;
rng(1342);

%% Instantiate the optimization problem.
N = 50;
theta = randn(3, N);
constants = [7; 1; 4; 8; 8];
obj = Example_1_Objective(constants);
cons = cell(N, 1);
for k = 1:N
    cons{k} = Example_1_Constraint(theta(:, k));
end
opt = Reduced_Space_Optimization_Under_Uncertainty(obj, cons);

%% Run finite difference checks.
opt.verbose = false;
z0 = rand(2, 1) + 1;
opt.Finite_Difference_Gradient_Check(z0);
opt.Finite_Difference_Hessian_Check(z0);

%% Do the optimization.
[u, z] = opt.Optimize(z0);

% save('reference_solution.mat','u','z')

u_ref = load('reference_solution.mat', 'u').u;
z_ref = load('reference_solution.mat', 'z').z;

error = max(norm(u - u_ref), norm(z - z_ref));
if error > 1e-12
    fprintf(2, '\noptimization_under_uncertainty/Example_1 failed.\n');
else
    fprintf(1, '\noptimization_under_uncertainty/Example_1 passed.\n');
end
