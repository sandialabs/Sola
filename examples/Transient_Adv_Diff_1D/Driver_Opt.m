clear
close all
clc
addpath(genpath('../../src'))

%% Set up the optimization problem.
m = 200;
N = 51;
T = 1;
num_space_control_nodes = 10;
n = num_space_control_nodes*(N-1);
obj = Adv_Diff_Gaussian_Source_Objective(m,n,T,N,num_space_control_nodes);
con = Adv_Diff_Gaussian_Source_Constraint(m,n,T,N,num_space_control_nodes);
opt = Reduced_Space_Optimization(obj,con);

%% Solve the optimization problem.
z0 = rand(n,1);
[u,z] = opt.Optimize(z0);

%% Compare the optimal state to the target.
x = con.x;
t = con.t_mesh;

u_reshape = reshape(u,m,N);
figure,
for k = 1:N
    target = obj.Evaluate_Target(t(k),x);
    plot(x,u_reshape(:,k),'-',x,target,'--','LineWidth',3)
    legend({'State','Target'})
    ylim([0 .2])
    pause(.05)
end
