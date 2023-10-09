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
adv_diff = Adv_Diff_Gaussian_Source(m,n,T,N,num_space_control_nodes);

%% Solve the optimization problem.
z0 = rand(n,1);
[u,z] = adv_diff.Optimize(z0);

%% Compare the optimal state to the target.
x = adv_diff.x;
t = adv_diff.t_mesh;

u_reshape = reshape(u,m,N);
figure,
for k = 1:N
    target = adv_diff.Evaluate_Target(t(k),x);
    plot(x,u_reshape(:,k),'-',x,target,'--','LineWidth',3)
    legend({'State','Target'})
    ylim([0 .2])
    pause(.05)
end
