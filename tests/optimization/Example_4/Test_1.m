%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
close all;
rng(89234);

n_y = 3;
n_z = 3;
T = .05;
n_t = 10^2;

obj = Example_4_Objective(n_y, n_z, T, n_t);
con = Example_4_Constraint(n_y, n_z, T, n_t);
opt = Reduced_Space_Optimization(obj, con);
con.verbose = false;
opt.verbose = false;

z0 = rand(n_z, 1) + 1;
y0 = rand(n_y, 1) + 1;
t = rand;
con.f_Jacobian_Check(y0, z0, t);
con.f_Hessian_Check(y0, z0, t);
opt.Finite_Difference_Gradient_Check(z0);
opt.Finite_Difference_Hessian_Check(z0);

[u, z] = opt.Optimize(z0);

% The optimal solution should be
% u(1:3:end) \approx exp(t)
% u(2:3:end) \approx exp(2*t)
% u(3:3:end) \approx exp(3*t)
% z \approx [1 ; 1 ; 1];

%%
u_sol = load('Solution_Example_4.mat', 'u').u;
z_sol = load('Solution_Example_4.mat', 'z').z;

error = 0;
error = max(error, norm(u_sol - u));
error = max(error, norm(z_sol - z));

if error > 1e-12
    fprintf(2, '\noptimization/Example_4 failed.\n');
else
    fprintf(1, '\noptimization/Example_4 passed.\n');
end

% save('Solution_Example_4.mat','u','z','obj')
