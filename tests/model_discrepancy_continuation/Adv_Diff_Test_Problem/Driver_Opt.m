clear
close all
clc
addpath(genpath('../../../src'))

m = 200;
diff_coeff = 1;
vel_coeff = 5;
robin_coeff = 2; 
reg_coeff = 10;
obj = Adv_Diff_Objective(m,reg_coeff);
con_hifi = Adv_Diff_Constraint(m,diff_coeff,vel_coeff,robin_coeff);
con_lofi = Diff_Constraint(obj,con_hifi);
opt_hifi = Reduced_Space_Optimization(obj,con_hifi);
opt_lofi = Reduced_Space_Optimization(obj,con_lofi);
x = con_hifi.x;

z0 = rand(m,1);
[u_hifi,z_hifi] = opt_hifi.Optimize(z0);
[u_lofi,z_lofi] = opt_lofi.Optimize(z0);

T = obj.T;
figure,
hold on
plot(x,u_hifi,'LineWidth',3)
plot(x,u_lofi,'--','LineWidth',3)
plot(x,T,'LineWidth',3)
xlabel('$x$','Interpreter','latex')
ylabel('Temperature','Interpreter','latex')
legend({'$u$','$\tilde{u}$','Target'},'location','south','Interpreter','latex')
set(gca, 'FontSize', 24); set(gcf, 'Color', 'White');

figure,
hold on
plot(x,z_hifi,'LineWidth',3)
plot(x,z_lofi,'--','LineWidth',3)
xlabel('$x$','Interpreter','latex')
ylabel('Source','Interpreter','latex')
legend({'$z$','$\tilde{z}$'},'location','north','Interpreter','latex')
set(gca, 'FontSize', 24); set(gcf, 'Color', 'White');

u_opt = u_lofi;
z_opt = z_lofi;

N = 5;
Z = zeros(m,N);
Z(:,1) = z_lofi;
E = (10^-2)*con_lofi.S + con_lofi.M;
for k = 2:N
    tmp = linsolve(E,randn(m,1)).^2;
    Z(:,k) = mean(z_lofi)*tmp/mean(tmp);
end
D = con_hifi.State_Solve(Z) - con_lofi.State_Solve(Z);

save('u_opt.mat','u_opt')
save('z_opt.mat','z_opt')
save('z_hifi.mat','z_hifi')
save('Z.mat','Z')
save('D.mat','D')
