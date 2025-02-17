%% Set up
clear;
close all;
clc;
addpath(genpath('../../src'));
load Optimization_Results.mat;

con_hifi = Adv_Diff_Constraint(m, diff_coeff, vel_coeff, robin_coeff);
con_lofi = Diff_Constraint(con_hifi);
x = con_lofi.x;

num_samples = 4;

E_z = (5.e-2) * con_hifi.S + con_hifi.M;

z = linsolve(E_z,sqrtm(con_hifi.M)*randn(m,num_samples-1));
z_lofi_norm = sqrt(z_lofi'*con_hifi.M*z_lofi);
z_norm = sqrt(diag(z'*con_hifi.M*z));
for k = 1:(num_samples-1)
    z(:,k) = z_lofi + .3*z_lofi_norm*z(:,k)/z_norm(k);
end
Z = zeros(m,num_samples);
Z(:,1) = z_lofi;
Z(:,2:end) = z;
D = con_hifi.State_Solve(Z) - con_lofi.State_Solve(Z);

data_interface = MD_Data_Interface_Discrepancy_Calibration(z_lofi,u_lofi,Z,D);
data_interface.Load_Data();

hyperparams = MD_Hyperparameters_Discrepancy_Calibration(data_interface,x);

u_prior_interface = MD_Analytic_Laplacian_u_Prior_Interface(con_hifi.M,hyperparams);
z_prior_interface = MD_Numeric_Laplacian_z_Prior_Interface(con_hifi.S,con_hifi.M,hyperparams);

alpha_u = (1.0) * hyperparams.alpha_u;
u_prior_interface.Set_alpha_u(alpha_u);
beta_u = (1.0) * hyperparams.beta_u;
u_prior_interface.Set_beta_u(beta_u);

alpha_z = (1.0) * hyperparams.alpha_z;
z_prior_interface.Set_alpha_z(alpha_z);
beta_z = (1.0) * hyperparams.beta_z;
z_prior_interface.Set_beta_z(beta_z);

hyperparams.Determine_alpha_d();
alpha_d = (1.0) * hyperparams.alpha_d;

%%
N = size(data_interface.Z,2);
colors = lines(N);

figure,
hold on
for k = 1:N
    plot(x,data_interface.Z(:,k),'LineWidth',3,'color',colors(k,:))
end

figure,
hold on
for k = 1:N
    plot(x,data_interface.D(:,k),'LineWidth',3,'color',colors(k,:))
end

%%
num_prior_samples = 100;
md_prior_sampling = MD_Prior_Sampling(data_interface, u_prior_interface, z_prior_interface);

delta_samples_at_zopt = md_prior_sampling.Prior_Discrepancy_Samples_at_z_opt(num_prior_samples);

figure;
hold on
plot(x, delta_samples_at_zopt, 'LineWidth', 3, 'color', [.9, .9, .9]);
plot(x, delta_samples_at_zopt(:, 1:10), 'LineWidth', 3);
ylim([-45,45])


%%
z = zeros(m, 2);
z(:, 1) = Z(:, 1) + .1 * x .* (1 - x);
z(:, 2) = Z(:, 1) .* (1 + .02 * cos(20 * pi * x));

figure;
hold on
plot(x, z(:,1), 'LineWidth', 3,'Color','magenta');
plot(x, z(:,2), 'LineWidth', 3,'Color','cyan');

[delta_samples,delta_zopt_samples] = md_prior_sampling.Prior_Discrepancy_Samples(z, num_prior_samples);
figure,
hold on
for k = 1:100
    plot(x, delta_samples{k}(:,1)-delta_zopt_samples(:,k), 'LineWidth', 3,'Color','magenta');
    plot(x, delta_samples{k}(:,2)-delta_zopt_samples(:,k), 'LineWidth', 3,'Color','cyan');
end
ylim([-45,45])



% %%
% md_post_samping = MD_Posterior_Sampling(data_interface, u_prior_interface, z_prior_interface);
% 
% num_post_samples = 500;
% md_post_samping.Compute_Posterior_Data(alpha_d, num_post_samples);
% Z_test = zeros(m, 3);
% Z_test(:, 1:2) = data_interface.Z;
% Z_test(:, 3) = 1.3 + 0.3*cos(4*pi*x);
% [delta_mean, delta_samples_at_zopt] = md_post_samping.Posterior_Discrepancy_Samples(Z_test);
% 
% figure,
% plot(x,delta_mean{1},x,D(:,1),'--','LineWidth',3)
% ylim([-20,5])
% legend({'$\delta(z_1,\overline{\theta})$','$d_1$'},'Position',[0.1768    0.7845    0.3161    0.1060],'Interpreter','latex')
% set(gca, 'fontsize', 18);
% 
% figure,
% plot(x,delta_mean{2},x,D(:,2),'--','LineWidth',3)
% ylim([-20,5])
% legend({'$\delta(z_2,\overline{\theta})$','$d_2$'},'Position',[0.1768    0.7845    0.3161    0.1060],'Interpreter','latex')
% set(gca, 'fontsize', 18);
% 
% figure,
% hold on
% plot(x,delta_samples_at_zopt{3},'color',[.8,.8,.8],'LineWidth',3)
% plot(x,delta_samples_at_zopt{3}(:,1:5),'LineWidth',3)
% ylim([-20,5])
% set(gca, 'fontsize', 18);