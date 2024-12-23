clear;
close all;
clc;
addpath(genpath('../../src'));
rng(1423);

plot_figures = true;
fd_check = true;

T = 30;
n_t = 100;
con = React_Rate_Eqn(T,n_t);
obj = Chem_React_Network_Objective(T, n_t,con);
opt = Reduced_Space_Optimization(obj, con);
SSA = SSA_System(con);

z0 = 2*rand;

if fd_check
    opt.Finite_Difference_Gradient_Check(z0);
    opt.Finite_Difference_Hessian_Check(z0);
end

[u_opt, z_opt] = opt.Optimize(z0);

num_samples = 500;
u_hifi = SSA.SSA_Mean(z_opt,num_samples);

y_opt = reshape(u_opt,9,n_t)'*con.nA*con.vol/con.state_scale;
y_hifi = reshape(u_hifi,9,n_t)'*con.nA*con.vol/con.state_scale;
t = linspace(0,T,n_t)';

if plot_figures
    for k = 1:9
        figure,
        plot(t,y_opt(:,k),t,y_hifi(:,k),'LineWidth',3)
        legend({'RRE','SSA'})
        title(['Species ',num2str(k)])
        set(gca,'FontWeight','Bold','FontSize',18)
    end
end

num_evals = 31;
z_range = z_opt*linspace(.95,1.05,num_evals)';
J_hifi = zeros(num_evals,1);
J_lofi = zeros(num_evals,1);
for k = 1:num_evals
    u_k = SSA.SSA_Mean(z_range(k),num_samples);
    J_hifi(k) = obj.J(u_k,z_range(k));
    J_lofi(k) = opt.Jhat(z_range(k));
end

figure,
plot(z_range,J_lofi,z_range,J_hifi,'LineWidth',3)
legend({'RRE','SSA'})
title('Objective Function')
set(gca,'FontWeight','Bold','FontSize',18)

