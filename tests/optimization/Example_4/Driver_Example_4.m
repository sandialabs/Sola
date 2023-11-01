clear;
close all;
addpath('../../../src/optimization/');
rng(89234);

m = 3;
n = 3;
T = .05;
N = 10^2;
obj = Example_4_Objective(m, n, T, N);
con = Example_4_Constraint(m, n, T, N);
opt = Reduced_Space_Optimization(obj, con);
con.verbose = false;
opt.verbose = false;
z0 = rand(n, 1) + 1;

y0 = rand(m, 1) + 1;
t = rand;
con.Time_Instance_RHS_Jacobian_y_Check(y0, z0, t);
con.Time_Instance_RHS_Jacobian_z_Check(y0, z0, t);
con.Time_Instance_RHS_Hessian_yy_Check(y0, z0, t);
con.Time_Instance_RHS_Hessian_yz_Check(y0, z0, t);
con.Time_Instance_RHS_Hessian_zy_Check(y0, z0, t);
con.Time_Instance_RHS_Hessian_zz_Check(y0, z0, t);
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
if error ~= 0
    disp('Error in example 4');
end

% save('Solution_Example_4.mat','u','z','obj')
