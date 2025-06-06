%% Set up
clear;
close all;
clc;
run('../../src/Set_Paths');

% Load data from surrogate optimization.
load Optimization_Results.mat;

% Set up high- and low-fidelity optimization problems.
obj = Adv_Diff_Objective(m, reg_coeff);
con_hifi = Adv_Diff_Constraint(m, diff_coeff, vel_coeff, robin_coeff);
con_lofi = Diff_Constraint(con_hifi);
opt_hifi = Reduced_Space_Optimization(obj, con_hifi);
opt_lofi = Reduced_Space_Optimization(obj, con_lofi);
x = con_lofi.x;

% Set up data and prior interfaces.
data_interface = MD_Data_Interface_Diff();
data_interface.Load_Data();

alpha_u = (1 / 2)^2;
alpha_z = (1 / 100)^2;
u_prior_interface = MD_Elliptic_u_Prior_Interface_Diff(alpha_u, opt_lofi);
z_prior_interface = MD_Elliptic_z_Prior_Interface_Diff(alpha_z, opt_lofi);

%%
md_post_sampling = MD_Posterior_Sampling(data_interface, u_prior_interface, z_prior_interface);
alpha_d = 1.e-2;
md_post_sampling.Compute_Posterior_Data(alpha_d, 1);
Z_test = Z;
[delta_mean, delta_samples] = md_post_sampling.Posterior_Discrepancy_Samples(Z_test);

%%
opt_prob_interface = MD_Opt_Prob_Interface_Sabl(opt_lofi, data_interface);
md_hessian_analysis = MD_Hessian_Analysis(opt_prob_interface, z_prior_interface);
md_update = MD_Update(md_post_sampling, md_hessian_analysis);

[z_update_mean, z_update_samples] = md_update.Posterior_Update_Samples();

%% Plot results.
figure;
hold on;
plot(x, z_lofi, 'color', 'black', 'LineWidth', 3);
plot(x, z_hifi, 'color', 'cyan', 'LineWidth', 3);
plot(x, z_update_mean, '--', 'color', 'red', 'LineWidth', 3);
legend({'Low-fidelity control', 'High-fidelity control', 'Update'});

u_hifi = con_hifi.State_Solve(z_hifi);
u_lofi = con_lofi.State_Solve(z_lofi);
u_before = con_hifi.State_Solve(z_lofi);
u_after = con_hifi.State_Solve(z_update_mean);
T = obj.T;

figure;
hold on;
plot(x, u_hifi, 'Color', 'cyan', 'LineWidth', 3);
plot(x, u_lofi, 'Color', 'black', 'LineWidth', 3);
plot(x, u_before, 'Color', 'magenta', 'LineWidth', 3);
plot(x, u_after, '--', 'Color', 'red', 'LineWidth', 3);
legend({'True model solution', 'Surrogate model solution', ...
        'True model with surrogate control', 'True model with updated control'});

% Export data for plotting externally.
save('plot_data.mat', 'x', 'z_hifi', 'z_lofi', 'z_update_mean', 'u_hifi', 'u_lofi', 'u_before', 'u_after', 'T');
