clear
close all
clc
addpath(genpath('../../src'))
rng(1423)

T = 10;
N = 100;
obj_hifi = Mass_Spring_Coupled(T,N);
obj_lofi = Mass_Spring_LoFi(obj_hifi);

z0 = 100*randn(N-1,1);

fd_check = false;
if fd_check
    obj_hifi.Finite_Difference_Gradient_Check(z0);
    obj_hifi.Finite_Difference_Hessian_Check(z0);
    
    obj_lofi.Finite_Difference_Gradient_Check(z0);
    obj_lofi.Finite_Difference_Hessian_Check(z0);
end

[u_lofi,z_lofi] = obj_lofi.Optimize(z0);
[~,z_hifi] = obj_hifi.Optimize(z0);

u_hifi = obj_hifi.State_Solve(z_lofi);
u_tmp = reshape(u_hifi,4,N)';
x1_hifi = u_tmp(:,1);
v1_hifi = u_tmp(:,2);
x2_hifi = u_tmp(:,3);
v2_hifi = u_tmp(:,4);
u_tmp = reshape(u_lofi,2,N)';
x1_lofi = u_tmp(:,1);
v1_lofi = u_tmp(:,2);
t = obj_hifi.t_mesh;

figure,
plot(t,x1_lofi,t,x1_hifi,'LineWidth',3)
xlabel('Time')
ylabel('$x_1$','Interpreter','latex')
legend({'Low-fidelity','High-fidelity'},'location','northwest')
set(gca,'fontsize', 18)

figure,
plot(t,v1_lofi,t,v1_hifi,'LineWidth',3)
xlabel('Time')
ylabel('$v_1$','Interpreter','latex')
legend({'Low-fidelity','High-fidelity'},'location','northwest')
set(gca,'fontsize', 18)

figure,
plot(t,0*t,t,x2_hifi,'LineWidth',3)
xlabel('Time')
ylabel('$x_2$','Interpreter','latex')
legend({'Low-fidelity','High-fidelity'},'location','northwest')
set(gca,'fontsize', 18)

figure,
plot(t,0*t,t,v2_hifi,'LineWidth',3)
xlabel('Time')
ylabel('$v_2$','Interpreter','latex')
legend({'Low-fidelity','High-fidelity'},'location','northwest')
set(gca,'fontsize', 18)

figure,
plot(t,obj_hifi.P_z*z_lofi,'LineWidth',3)
xlabel('Time')
ylabel('$z$','Interpreter','latex')
set(gca,'fontsize', 18)

Z = zeros(N-1,2);
Z(:,1) = z_lofi;
Z(:,2) = 100*t(2:end);

D = zeros(2*N,2);
for k = 1:size(D,2)
    tmp = obj_hifi.State_Solve(Z(:,k));
    tmp = reshape(tmp,4,N)';
    tmp = tmp(:,1:2)';
    D(:,k) = tmp(:) - obj_lofi.State_Solve(Z(:,k));
end

save('Optimization_Results.mat','z_lofi','z_hifi','u_lofi','Z','D','T')