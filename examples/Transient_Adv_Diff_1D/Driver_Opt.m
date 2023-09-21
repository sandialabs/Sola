clear
close all
clc
addpath(genpath('../../src'))

m = 200;
diff_coeff = 1;
vel_coeff = 1/2;
robin_coeff = 3/2; 
reg_coeff = 10;
obj_hifi = Adv_Diff(m,diff_coeff,vel_coeff,robin_coeff,reg_coeff);
obj_lofi = Diff(obj_hifi);

z0 = rand(m,1);
[u_hifi,z_hifi] = obj_hifi.Optimize(z0);
[u_lofi,z_lofi] = obj_lofi.Optimize(z0);

x = obj_hifi.x;
T = obj_hifi.T;
u = obj_hifi.State_Solve(z_lofi);
ymax = 1.1*max([u_lofi;u_hifi;u;T]);

figure,
hold on
plot(x,T,'LineWidth',3)
plot(x,u,'LineWidth',3)
plot(x,u_lofi,'LineWidth',3)
xlabel('$x$','Interpreter','latex')
ylabel('Concentration','Interpreter','latex')
ylim([0,ymax])
title('For low-fidelity control','Interpreter','latex')
legend({'Target','$u$','$\tilde{u}$'},'location','south','Interpreter','latex')
set(gca, 'FontSize', 24); set(gcf, 'Color', 'White');

figure,
hold on
plot(x,T,'LineWidth',3)
plot(x,u_hifi,'LineWidth',3)
xlabel('$x$','Interpreter','latex')
ylabel('Concentration','Interpreter','latex')
ylim([0,ymax])
title('For high-fidelity control','Interpreter','latex')
legend({'Target','$u$'},'location','south','Interpreter','latex')
set(gca, 'FontSize', 24); set(gcf, 'Color', 'White');

figure,
hold on
plot(x,z_hifi,'LineWidth',3,'color',[0.8500 0.3250 0.0980])
plot(x,z_lofi,'LineWidth',3,'color',[0.9290 0.6940 0.1250])
xlabel('$x$','Interpreter','latex')
ylabel('Source','Interpreter','latex')
legend({'$z$','$\tilde{z}$'},'location','north','Interpreter','latex')
set(gca, 'FontSize', 24); set(gcf, 'Color', 'White');

Z = zeros(m,2);
Z(:,1) = z_lofi;
Z(:,2) = ones(m,1);

Y = zeros(m,2);
for k = 1:size(Y,2)
    Y(:,k) = obj_hifi.State_Solve(Z(:,k)) - obj_lofi.State_Solve(Z(:,k));
end

save('Optimization_Results.mat','m','diff_coeff','vel_coeff','robin_coeff','reg_coeff','z_lofi','z_hifi','u_lofi','Z','Y')