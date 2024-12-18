clear;
close all;
addpath('../../../src/optimization/');
rng(132);

%% Set up problem.

n_y = 2;
T = .1;
n_t = 5 * 10^2;
n_z = n_t - 1;

obj = Example_5_Objective(n_y, n_z, T, n_t);
con = Example_5_Constraint(n_y, n_z, T, n_t);
opt = Reduced_Space_Optimization(obj, con);
opt.verbose = false;

%% Finite difference checks.

z0 = rand(n_z, 1);
% u0 = rand(n_y * n_t, 1);
% obj.Finite_Difference_Gradient_Check(u0, z0);
% obj.Finite_Difference_Hessian_Check(u0, z0);

[u, z] = opt.Optimize(z0);

%% Check against saved solution.

% % The optimal solution should be
% % u(1:2:end) \approx exp(t^2)
% % u(2:2:end) \approx exp(t^3)
% % z \approx obj.z_time_mesh.^2

u_sol = load('Solution_Example_5.mat', 'u').u;
z_sol = load('Solution_Example_5.mat', 'z').z;

error = 0;
error = max(error, norm(u_sol - u)/norm(u_sol));
error = max(error, norm(z_sol - z)/norm(z_sol));
if error > 1.e-6
    disp('Error in example 5');
    disp('Computed objective:');
    disp(num2str(obj.J(u, z), '%.8e'));
    disp('Saved objective:');
    disp(num2str(obj.J(u_sol, z_sol), '%.8e'));
end

% save('Solution_Example_5.mat','u','z','obj')
