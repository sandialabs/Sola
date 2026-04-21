%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
close all;
rng(142534);

n_y = 2;
n_z = 2;
T = .1;
n_t = 10^2;

obj = Example_2_Objective(n_y, n_z, T, n_t);
con = Example_2_Constraint(n_y, n_z, T, n_t);
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

if error > 1e-12
    fprintf(2, '\noptimization/Example_2 failed.\n');
else
    fprintf(1, '\noptimization/Example_2 passed.\n');
end

% save('Solution_Example_2.mat','u','z','obj')
