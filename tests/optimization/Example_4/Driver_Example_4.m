clear
close all
addpath('../../../src/Optimization/')
rng(89234)

m = 3;
n = 3;
T = .05;
N = 10^2;
obj = Example_4(m,n,T,N);
obj.verbose = false;
z0 = rand(n,1)+1;

u0 = rand(m,1)+1;
t = rand;
obj.Time_Instance_RHS_Jacobian_u_Check(u0,z0,t);
obj.Time_Instance_RHS_Jacobian_z_Check(u0,z0,t);
obj.Time_Instance_RHS_Hessian_uu_Check(u0,z0,t);
obj.Time_Instance_RHS_Hessian_uz_Check(u0,z0,t);
obj.Time_Instance_RHS_Hessian_zu_Check(u0,z0,t);
obj.Time_Instance_RHS_Hessian_zz_Check(u0,z0,t);
obj.Finite_Difference_Gradient_Check(z0);
obj.Finite_Difference_Hessian_Check(z0);

[u,z] = obj.Optimize(z0);

% The optimal solution should be
% u(1:3:end) \approx exp(t)
% u(2:3:end) \approx exp(2*t)
% u(3:3:end) \approx exp(3*t)
% z \approx [1 ; 1 ; 1];

%%
u_sol = load('Solution_Example_4.mat','u').u;
z_sol = load('Solution_Example_4.mat','z').z;

error = 0;
error = max(error,norm(u_sol-u));
error = max(error,norm(z_sol-z));
if error ~= 0
   disp('Error in example 4') 
end

% save('Solution_Example_4.mat','u','z','obj')