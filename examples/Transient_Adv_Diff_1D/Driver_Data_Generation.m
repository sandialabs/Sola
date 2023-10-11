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
con = Adv_Diff_Gaussian_Source_Constraint(m,n,T,N,num_space_control_nodes);

%% Solve the state equation for several random controls.
num_solves = 3;
Z = randn(n,num_solves);
Y = zeros(m*N,num_solves);
for k = 1:num_solves
    Y(:,k) = con.State_Solve(Z(:,k));
end

%% Visualize each set of drawn states.
x = con.x;
t = con.t_mesh;
for j = 1:num_solves
    u_reshape = reshape(Y(:,j),m,N);
    figure,
    for k = 1:N
        plot(x,u_reshape(:,k),'LineWidth',3)
        ylim([-.1 .1])
        title(['Forward solve ',num2str(j),' at time ',num2str(t(k))])
        pause(.05)
    end
end
