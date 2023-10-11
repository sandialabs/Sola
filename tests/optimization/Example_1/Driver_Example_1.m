%%
% Clear workspace and add path
clear
close all
addpath('../../../src/Optimization/')
rng(1342)

%%
% Instantiate the Example_1 object
obj = Example_1_Objective();
con = Example_1_Constraint();
opt = Reduced_Space_Optimization(obj,con);
opt.verbose = false;

%%
% Generate a random control and execute finite difference tests
z0 = rand(2,1)+1;
opt.Finite_Difference_Gradient_Check(z0);
opt.Finite_Difference_Hessian_Check(z0);

%%
% Execute optimization
% The optimal solution should be
% u = ( 7  1  4)
% z = (8  8)
[u,z] = opt.Optimize(z0);

%%
u_sol = load('Solution_Example_1.mat','u').u;
z_sol = load('Solution_Example_1.mat','z').z;

error = 0;
error = max(error,norm(u_sol-u));
error = max(error,norm(z_sol-z));
if error ~= 0
   disp('Error in example 1') 
end

% save('Solution_Example_1.mat','u','z','obj')