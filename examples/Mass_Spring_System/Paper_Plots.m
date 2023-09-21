clear
close all
clc

addpath(genpath('../../src'))
load('HDSA_Results.mat')

working_path = pwd;
write_path = '~/Desktop/Model_Discrepancy_Sampling/figures/mass_spring/';

figure,
hold on
plot(t(2:end),z_lofi,'LineWidth',3,'color','black')
xlabel('$t$','Interpreter','latex')
ylabel('Forcing')
legend('$\tilde{z}$','Interpreter','latex','FontSize',18,'Location','northwest')
set(gca,'fontsize', 18)
cd(write_path)
saveas(gca,'lofi_forcing','epsc')
cd(working_path)

T = obj_hifi.target(t);
u_low = obj_lofi.State_Solve(z_lofi);
u_high = obj_hifi.State_Solve(z_lofi);
figure,
hold on
plot(t,T,':','LineWidth',3,'color','magenta')
plot(t,u_low(1:2:end),'LineWidth',3,'color','green')
plot(t,u_high(1:4:end),'LineWidth',3,'color','black')
plot(t,T,'--','LineWidth',3,'color','magenta')
xlabel('$t$','Interpreter','latex')
ylabel('Block 1 displacement')
legend({'$T$','$\tilde{S}(\tilde{z})$','$S(\tilde{z})$'},'Interpreter','latex','FontSize',18,'Location','northwest')
set(gca,'fontsize', 18)
cd(write_path)
saveas(gca,'states_at_lofi_forcing','epsc')
cd(working_path)

u_true_update = obj_hifi.State_Solve(z_update_mean);
figure,
hold on
plot(t(2:end),z_lofi,'color','black','LineWidth',3)
plot(t(2:end),z_hifi,'color','cyan','LineWidth',3)
plot(t(2:end),z_update_mean,'--','color','red','LineWidth',3)
for k = 1:num_post_samples
    plot(t(2:end),z_update_samples(:,k),'color',[.9,.9,.9],'LineWidth',3)
end
plot(t(2:end),md_update.z_opt,'color','black','LineWidth',3)
plot(t(2:end),z_hifi,'color','cyan','LineWidth',3)
plot(t(2:end),z_update_mean,'--','color','red','LineWidth',3)
xlabel('$t$','Interpreter','latex')
ylabel('Forcing')
legend({'$\tilde{z}$','$z^\star$','$\overline{z}$'},'Location','northwest','Interpreter','latex')
set(gca,'fontsize', 18)
cd(write_path)
saveas(gca,'post_opt_solution','epsc')
cd(working_path)

u_true_update_samples = zeros(length(u_high),num_post_samples);
for k = 1:num_post_samples
   u_true_update_samples(:,k) = obj_hifi.State_Solve(z_update_samples(:,k)); 
end

u_high_high = obj_hifi.State_Solve(z_hifi);
figure,
hold on
plot(t,u_high(1:4:end),'LineWidth',3,'color','black')
plot(t,u_high_high(1:4:end),'LineWidth',3,'color','cyan')
plot(t,u_true_update(1:4:end),'--','LineWidth',3,'color','red')
for k = 1:num_post_samples
    plot(t,u_true_update_samples(1:4:end,k),'color',[.9,.9,.9],'LineWidth',3)
end
plot(t,u_high(1:4:end),'LineWidth',3,'color','black')
plot(t,u_high_high(1:4:end),'LineWidth',3,'color','cyan')
plot(t,u_true_update(1:4:end),'--','LineWidth',3,'color','red')
xlabel('$t$','Interpreter','latex')
ylabel('Block 1 displacement')
legend({'$S(\tilde{z})$','$S(z^\star)$','$S(\overline{z})$'},'Interpreter','latex','FontSize',18,'Location','northwest')
set(gca,'fontsize', 18)
cd(write_path)
saveas(gca,'post_opt_solution_states','epsc')
cd(working_path)