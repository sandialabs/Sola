clear
close all
clc
addpath(genpath('../../../src'))

m = 200;
diff_coeff = 1;
vel_coeff = 1/2;
robin_coeff = 2; 
reg_coeff = 10;
obj_hifi = Adv_Diff(m,diff_coeff,vel_coeff,robin_coeff,reg_coeff);
obj_lofi = Diff(obj_hifi);
x = obj_hifi.x;

z0 = rand(m,1);
[u_hifi,z_hifi] = obj_hifi.Optimize(z0);
[u_lofi,z_lofi] = obj_lofi.Optimize(z0);

T = obj_hifi.T;
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
