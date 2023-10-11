clear
close all
addpath('../../../src/Optimization/')
rng(132)

m = 2;
T = .1;
N = 5*10^2;
n = N-1;
obj = Example_5_Objective(m,n,T,N);
con = Example_5_Constraint(m,n,T,N);
opt = Reduced_Space_Optimization(obj,con);
opt.verbose = false;
z0 = rand(n,1);
%obj.Finite_Difference_Gradient_Check(z0);
%obj.Finite_Difference_Hessian_Check(z0);
[u,z] = opt.Optimize(z0);

% % The optimal solution should be
% % u(1:2:end) \approx exp(t^2)
% % u(2:2:end) \approx exp(t^3)
% % z \approx obj.z_time_mesh.^2
%
%%
u_sol = load('Solution_Example_5.mat','u').u;
z_sol = load('Solution_Example_5.mat','z').z;

error = 0;
error = max(error,norm(u_sol-u));
error = max(error,norm(z_sol-z));
if error ~= 0
   disp('Error in example 5')
   disp('Computed objective:')
   disp(num2str(opt.Objective(u,z), '%.8e'))
   disp('Saved objective:')
   disp(num2str(opt.Objective(u_sol,z_sol), '%.8e'))
end

% save('Solution_Example_5.mat','u','z','obj')
