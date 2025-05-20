clear;
close all;
addpath(genpath('../../../src'));
rng(132253);

random_numbers = randn(10^6, 1);
writematrix(random_numbers, 'random_numbers.txt');

rng(132253);

n_y = 50;
n_t = 11;
T = 1;
n_z = n_y;
obj = Adv_Diff_Objective(n_y, n_z, T, n_t);
obj.w(1) = 0;
obj.w(end) = obj.w(end - 1);
con_hifi = Adv_Diff_Constraint(n_y, n_z, T, n_t);
con = Diff_Constraint(n_y, n_z, T, n_t);
opt = Reduced_Space_Optimization(obj, con);
x = con.x;
t = linspace(0, T, n_t)';

%%
data_interface = MD_Data_Interface_Transient_Test_Problem();
data_interface.Load_Data();
opt_prob_interface = MD_Opt_Prob_Interface_Sabl(opt, data_interface);

writematrix(data_interface.u_opt, 'u_opt.txt');
writematrix(data_interface.z_opt, 'z_opt.txt');
writematrix(data_interface.D, 'D.txt');
writematrix(data_interface.Z, 'Z.txt');

adapt_time_variance = true;
u_hyperparam_interface = MD_u_Hyperparameter_Interface_Transient_Test_Problem(x, t, adapt_time_variance);
u_hyperparam_interface.alpha_u = 1;
u_hyperparam_interface.beta_u = .009;
u_hyperparam_interface.beta_t = .026;
u_hyperparam_interface.alpha_d = 2e-6;
u_hyperparam_interface.gsvd_num_sing_vals = 50;
u_hyperparam_interface.gsvd_oversampling = 0;
u_hyperparam_interface.gsvd_num_subspace_iter = 1;

spatial_u_prior_interface = MD_Numeric_Laplacian_u_Prior_Interface(con.S, con.M, data_interface, u_hyperparam_interface);
transient_prior_cov = MD_Transient_Prior_Covariance_Sabl(data_interface, u_hyperparam_interface, T, n_t, n_y);
u_prior_interface = MD_Transient_Elliptic_u_Prior_Interface(data_interface, spatial_u_prior_interface, transient_prior_cov);

z_hyperparam_interface = MD_z_Hyperparameter_Interface_Transient_Test_Problem(x, con);
z_hyperparam_interface.alpha_z = .7;
z_hyperparam_interface.beta_z = .004;
z_prior_interface = MD_Numeric_Laplacian_z_Prior_Interface(con.S, con.M, data_interface, z_hyperparam_interface, u_prior_interface);

%%
num_prior_samples = 100;
md_prior_sampling = MD_Prior_Sampling(data_interface, u_prior_interface, z_prior_interface);
md_prior_viz = MD_Prior_Visualization(md_prior_sampling);
md_prior_sampling.Generate_Prior_Discrepancy_z_opt_Sample_Data(num_prior_samples);
md_prior_sampling.Generate_Prior_Discrepancy_z_pert_Sample_Data();
%md_prior_viz.Visualization_for_Prior_Time_Evolution(1,false);
%md_prior_viz.Visualization_for_Prior_Discrepancy_at_z_opt(1);
%md_prior_viz.Visualization_for_Prior_Discrepancy_at_z_pert(1);

%%
md_post_sampling = MD_Posterior_Sampling(data_interface, u_prior_interface, z_prior_interface);
num_post_samples = 100;
md_post_sampling.Compute_Posterior_Data(u_hyperparam_interface.alpha_d, num_post_samples);

%%
md_hessian_analysis = MD_Hessian_Analysis(opt_prob_interface, z_prior_interface);
num_evals = 5;
oversampling = 4;
md_hessian_analysis.Compute_Hessian_GEVP(data_interface.z_opt, num_evals, oversampling);
md_update = MD_Update(md_post_sampling, md_hessian_analysis);

[z_update_mean, z_update_samples] = md_update.Posterior_Update_Samples();

%%
opt_hifi = Reduced_Space_Optimization(obj, con_hifi);
[u_hifi, z_hifi] = opt_hifi.Optimize(data_interface.z_opt);

figure;
hold on;
plot(x, data_interface.z_opt, 'color', 'black', 'LineWidth', 3);
plot(x, z_hifi, 'color', 'cyan', 'LineWidth', 3);
plot(x, z_update_mean, '--', 'color', 'red', 'LineWidth', 3);
for k = 1:num_post_samples
    plot(x, z_update_samples(:, k), 'color', [.9, .9, .9], 'LineWidth', 3);
end
plot(x, data_interface.z_opt, 'color', 'black', 'LineWidth', 3);
plot(x, z_hifi, 'color', 'cyan', 'LineWidth', 3);
plot(x, z_update_mean, '--', 'color', 'red', 'LineWidth', 3);
