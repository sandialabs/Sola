%% Set up
clear;
close all;
clc;
addpath(genpath('../../src'));
load Optimization_Results.mat;

con_hifi = Adv_Diff_Constraint(m, diff_coeff, vel_coeff, robin_coeff);
con_lofi = Diff_Constraint(con_hifi);
x = con_lofi.x;

data_interface = MD_Data_Interface_Discrepancy_Calibration();
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
md_post_samping = MD_Posterior_Sampling(data_interface, u_prior_interface, z_prior_interface);

num_post_samples = 500;
md_post_samping.Compute_Posterior_Data(alpha_d, num_post_samples);
Z_test = zeros(m, 3);
Z_test(:, 1:2) = data_interface.Z;
Z_test(:, 3) = 1.3 + 0.3*cos(4*pi*x);
[delta_mean, delta_samples] = md_post_samping.Posterior_Discrepancy_Samples(Z_test);

figure,
plot(x,delta_mean{1},x,D(:,1),'--','LineWidth',3)
ylim([-20,5])
legend({'$\delta(z_1,\overline{\theta})$','$d_1$'},'Position',[0.1768    0.7845    0.3161    0.1060],'Interpreter','latex')
set(gca, 'fontsize', 18);

figure,
plot(x,delta_mean{2},x,D(:,2),'--','LineWidth',3)
ylim([-20,5])
legend({'$\delta(z_2,\overline{\theta})$','$d_2$'},'Position',[0.1768    0.7845    0.3161    0.1060],'Interpreter','latex')
set(gca, 'fontsize', 18);

figure,
hold on
plot(x,delta_samples{3},'color',[.8,.8,.8],'LineWidth',3)
plot(x,delta_samples{3}(:,1:5),'LineWidth',3)
ylim([-20,5])
set(gca, 'fontsize', 18);