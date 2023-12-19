%%
clear;
close all;
clc;
run('../../src/Set_Paths');

con_hifi = load('HiFi_Opt_Results.mat', 'con').con;
reg_coeff = load('HiFi_Opt_Results.mat', 'obj').obj.reg_coeff;
con = Thermochemical_LoFi_Constraint_AD(con_hifi);
con.AD_Initialization('Lofi_AD_Files');
obj = Thermochemical_Dynamic_Objective(con, reg_coeff);
opt = Reduced_Space_Optimization(obj, con);

n_y = con.fe.m;
n_t = con.n_t;
T = con.T;
control_time_nodes = con.control_time_nodes;

%%
data_interface = Thermochemical_Data_Interface(n_y, n_t);
data_interface.Load_Data();

opt_prob_interface = MD_Opt_Prob_Interface_Sabl(opt, data_interface);

beta_tu = 10;
beta_iu = 1.e5;
beta_td = 0.0;
transient_prior_cov = MD_Transient_Prior_Covariance_Sabl(beta_tu, beta_iu, beta_td, T, n_t, 4 * n_y);

alpha_u = 200.0^2;
u_prior_interface = Thermochemical_Elliptic_u_Prior_Interface(alpha_u, transient_prior_cov, con.fe);

alpha_z = 1.e-10;
z_prior_interface = Thermochemical_Elliptic_z_Prior_Interface(alpha_z, con.fe, T, control_time_nodes);

%%
md_prior_sampling = MD_Prior_Sampling(data_interface, u_prior_interface, z_prior_interface);

num_prior_samples = 10;
prior_delta_samples = md_prior_sampling.Prior_Discrepancy_Samples_at_z_opt(num_prior_samples);

fig_nums = 1:4;
for k = 1:num_prior_samples
    Plot_States(prior_delta_samples(:, k), con, fig_nums);
    pause();
end

%%
md_hessian_analysis = MD_Hessian_Analysis(opt_prob_interface, z_prior_interface);
md_update = MD_Update(opt_prob_interface, data_interface, u_prior_interface, z_prior_interface, md_hessian_analysis);

alpha_d = 1.e-4;
num_post_samples = 5;
md_update.Compute_Posterior_Data(alpha_d, num_post_samples);

%%
Z_test = zeros(n_y * control_time_nodes, 2);
Z_test(:, 1) = data_interface.Z;
Z_test(:, 2) = data_interface.Z + 20;
[delta_mean, delta_samples] = md_update.Posterior_Discrepancy_Samples(Z_test);

Plot_States(data_interface.D, con, fig_nums);
fig_nums = 4 + fig_nums;
Plot_States(delta_mean{1}, con, fig_nums);
fig_nums = 4 + fig_nums;
Plot_States(delta_mean{2}, con, fig_nums);

%%
num_evals = 51;
oversampling = 20;
md_update.md_hessian_analysis.Compute_Hessian_GEVP(data_interface.z_opt, num_evals, oversampling);

figure;
plot(md_update.md_hessian_analysis.evals, 'o');
title('Hessian Eigenvalues');
set(gca, 'fontsize', 18);

[z_update_mean, z_update_samples] = md_update.Posterior_Update_Samples();

save('MD_Analysis.mat');
