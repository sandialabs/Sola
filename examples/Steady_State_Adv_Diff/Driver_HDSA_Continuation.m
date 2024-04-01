%% Set up
clear;
close all;
clc;
addpath(genpath('../../src'));
load Optimization_Results.mat;

suppress_figures = false;

obj = Adv_Diff_Objective(m, reg_coeff);
con_hifi = Adv_Diff_Constraint(m, diff_coeff, vel_coeff, robin_coeff);
con_lofi = Diff_Constraint(con_hifi);
opt_hifi = Reduced_Space_Optimization(obj, con_hifi);
opt_lofi = Reduced_Space_Optimization(obj, con_lofi);
x = con_lofi.x;

data_interface = MD_Data_Interface_Diff(3, 5, 'OED');
data_interface.Load_Data();
opt_prob_interface = MD_Opt_Prob_Interface_Sabl(opt_lofi, data_interface);
alpha_u = (1 / 2)^2;
alpha_z = (1 / 100)^2;
u_prior_interface = MD_Elliptic_u_Prior_Interface_Diff(alpha_u, opt_lofi);
z_prior_interface = MD_Elliptic_z_Prior_Interface_Diff(alpha_z, opt_lofi);

%%
md_hessian_analysis = MD_Hessian_Analysis(opt_prob_interface, z_prior_interface);
num_continuation_steps = 50;
md_update = MD_Continuation_Update(opt_prob_interface, data_interface, u_prior_interface, z_prior_interface, md_hessian_analysis, num_continuation_steps);

alpha_d = 1.e-2;
num_post_samples = 100;
md_update.Compute_Posterior_Data(alpha_d, num_post_samples);
[u_update_mean, z_update_mean] = md_update.Posterior_Update_Mean();

if ~suppress_figures
    figure;
    hold on;
    plot(x, z_lofi, 'color', 'black', 'LineWidth', 3);
    plot(x, z_hifi, 'color', 'cyan', 'LineWidth', 3);
    plot(x, z_update_mean(:, end), '--', 'color', 'red', 'LineWidth', 3);
    legend({'Low-fidelity control', 'High-fidelity control', 'Update'});
end
