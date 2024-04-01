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
md_oed.Offline_Computation();

samps_per_N = 5;
N_range = (2:4)';
p = length(N_range);

oed_beta_samps = cell(p, samps_per_N);
oed_Z_samps = cell(p, samps_per_N);
oed_D_samps = cell(p, samps_per_N);
oed_obj = zeros(p, samps_per_N);

rand_Z_samps = cell(p, samps_per_N);
rand_D_samps = cell(p, samps_per_N);
subrand_Z_samps = cell(p, samps_per_N);
subrand_D_samps = cell(p, samps_per_N);

for k = 1:p
    N = N_range(k);
    for i = 1:samps_per_N
        beta_0 = randn(num_evals * (N - 1), 1);
        alpha_d = 1.e-7;
        reg_coeff = 1.e-5;
        [oed_beta_samps{k, i}, oed_Z_samps{k, i}] = md_oed.Generate_Optimal_Design(beta_0, alpha_d, reg_coeff);
        oed_D_samps{k, i} = zeros(length(u_lofi), N);
        for j = 1:N
            oed_D_samps{k, i}(:, j) = con_hifi.State_Solve(oed_Z_samps{k, i}(:, j)) - con_lofi.State_Solve(oed_Z_samps{k, i}(:, j));
        end
        oed_obj(k, i) = md_oed.Evaluate_OED_Objective(oed_beta_samps{k, i}, alpha_d, reg_coeff);

        rand_Z_samps{k, i} = md_oed.Generate_Random_Design(N);
        rand_D_samps{k, i} = zeros(length(u_lofi), N);
        for j = 1:N
            rand_D_samps{k, i}(:, j) = con_hifi.State_Solve(rand_Z_samps{k, i}(:, j)) - con_lofi.State_Solve(rand_Z_samps{k, i}(:, j));
        end
        %         subrand_Z_samps{k, i} = md_oed.Generate_Random_Design_from_Subspace(N);
        %         subrand_D_samps{k, i}  = zeros(length(u_lofi),N);
        %         for j = 1:N
        %             subrand_D_samps{k, i}(:,j) = con_hifi.State_Solve(subrand_Z_samps{k, i}(:,j)) - con_lofi.State_Solve(subrand_Z_samps{k, i}(:,j));
        %         end
    end
end

save('OED_Ensemble_Results.mat');
