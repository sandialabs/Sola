clear
close all
clc
addpath(genpath('../../src'))

m = 200;
N = 51;
T = 1;
num_space_control_nodes = 10;
n = num_space_control_nodes*(N-1);
adv_diff = Adv_Diff_Gaussian_Source(m,n,T,N,num_space_control_nodes);

x = adv_diff.x;
t = adv_diff.t_mesh;
T = adv_diff.T;

num_solves = 3;
Z = randn(n,num_solves);
Y = zeros(m*N,num_solves);
for k = 1:num_solves
    Y(:,k) = adv_diff.State_Solve(Z(:,k));
end

for j = 1:num_solves
    u_reshape = reshape(Y(:,j),m,N);
    figure,
    for k = 1:N
        plot(x,u_reshape(:,k),'LineWidth',3)
        title(['Forward solve ',num2str(j),' at time ',num2str(t(k))])
        pause(.05)
    end
end
