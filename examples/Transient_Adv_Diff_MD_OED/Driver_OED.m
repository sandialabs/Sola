%% Set up
clear;
close all;
clc;
addpath(genpath('../../src'));
load Optimization_Results.mat;

suppress_figures = false;

n_z = num_space_control_nodes * (n_t - 1);
obj = Adv_Diff_Gaussian_Source_Objective(n_y, n_z, T, n_t, num_space_control_nodes, reg_coeff);
con_lofi = Adv_Diff_Gaussian_Source_Constraint(n_y, n_z, T, n_t, num_space_control_nodes, diff_coeff, vel_coeff_lofi);
opt_lofi = Reduced_Space_Optimization(obj, con_lofi);
con_hifi = Adv_Diff_Gaussian_Source_Constraint(n_y, n_z, T, n_t, num_space_control_nodes, diff_coeff, vel_coeff_hifi);
opt_hifi = Reduced_Space_Optimization(obj, con_hifi);
x = con_lofi.x;

data_interface = MD_Data_Interface_Adv_Diff();
data_interface.Load_Data();

beta_t = 50;
beta_i = 1.e5;
transient_prior_cov = MD_Transient_Prior_Covariance_Sabl(beta_t, beta_i, T, n_t, n_y);

alpha_u = (1 / 1)^2;
u_prior_interface = MD_Transient_Elliptic_u_Prior_Interface_Adv_Diff(alpha_u, transient_prior_cov, opt_lofi);

z_prior_interface = MD_z_Prior_Interface_Adv_Diff(obj);

%%
opt_prob_interface = MD_Opt_Prob_Interface_Sabl(opt_lofi, data_interface);
md_hessian_analysis = MD_Hessian_Analysis(opt_prob_interface, z_prior_interface);
num_evals = 26;
oversampling = 20;
md_hessian_analysis.Compute_Hessian_GEVP(data_interface.z_opt, num_evals, oversampling);

%% OED
oed_interface = MD_OED_Interface_Adv_Diff(data_interface, obj);

md_oed = MD_OED(opt_prob_interface, data_interface, u_prior_interface, z_prior_interface, md_hessian_analysis, oed_interface);

N = 5;

run_random = false;
if run_random
    Z = md_oed.Generate_Random_Design(N);
    for k = 2:N
        Plot_Control(x, Z(:, k), con_lofi);
    end
end

md_oed.Offline_Computation();

run_random_subspace = false;
if run_random_subspace
    Z = md_oed.Generate_Random_Design_from_Subspace(N);
    for k = 2:N
        Plot_Control(x, Z(:, k), con_lofi);
    end
end

run_oed = false;
if run_oed
    beta_0 = randn(num_evals * (N - 1), 1);
    alpha_d = 1.e-7;
    reg_coeffs = 10.^(-4:-1:-8)';
    [beta_L_curve, Z_L_curve, post_var, reg_val] = md_oed.L_Curve_Analysis(beta_0, alpha_d, reg_coeffs);
    figure;
    plot(post_var, reg_val, 'o', 'MarkerSize', 10);
    set(gca, 'fontsize', 18);

    L_curve_index = 2;
    Z = Z_L_curve{L_curve_index};

    for k = 2:N
        Plot_Control(x, Z(:, k), con_lofi);
    end
end
