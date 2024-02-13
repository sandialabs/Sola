%% Set up
clear;
close all;
clc;
addpath(genpath('../../src'));
load Optimization_Results.mat;

obj = Adv_Diff_Objective(m, reg_coeff);
con_hifi = Adv_Diff_Constraint(m, diff_coeff, vel_coeff, robin_coeff);
con_lofi = Diff_Constraint(con_hifi);
opt_hifi = Reduced_Space_Optimization(obj, con_hifi);
opt_lofi = Reduced_Space_Optimization(obj, con_lofi);
x = con_lofi.x;

%% HDSA interfaces
data_interface = MD_Data_Interface_Diff();
data_interface.Load_Data();
opt_prob_interface = MD_Opt_Prob_Interface_Sabl(opt_lofi, data_interface);
alpha_u = (1 / 2)^2;
alpha_z = (1 / 100)^2;
u_prior_interface = MD_Elliptic_u_Prior_Interface_Diff(alpha_u, opt_lofi);
z_prior_interface = MD_Elliptic_z_Prior_Interface_Diff(alpha_z, opt_lofi);

%% Hessian analysis
md_hessian_analysis = MD_Hessian_Analysis(opt_prob_interface, z_prior_interface);
num_evals = 20;
oversampling = 10;
md_hessian_analysis.Compute_Hessian_GEVP(data_interface.z_opt, num_evals, oversampling);

%% OED
oed_interface = MD_OED_Interface_Diff(data_interface, con_lofi);

md_oed = MD_OED(opt_prob_interface, data_interface, u_prior_interface, z_prior_interface, md_hessian_analysis, oed_interface);

md_oed.Offline_Computation();

samps_per_N = 20;
N_range = (2:6)';
p = length(N_range);

oed_beta_samps = cell(p, samps_per_N);
oed_Z_samps = cell(p, samps_per_N);
oed_D_samps = cell(p, samps_per_N);

rand_Z_samps = cell(p, samps_per_N);
rand_D_samps = cell(p, samps_per_N);
subrand_Z_samps = cell(p, samps_per_N);
subrand_D_samps = cell(p, samps_per_N);

for k = 1:p
    N = N_range(k);
    for i = 1:samps_per_N
        beta_0 = randn(num_evals * (N - 1), 1);
        alpha_d = 1.e-2;
        reg_coeff = 1.e-6;
        [oed_beta_samps{k, i}, oed_Z_samps{k, i}] = md_oed.Generate_Optimal_Design(beta_0, alpha_d, reg_coeff);
        oed_D_samps{k, i} = con_hifi.State_Solve(oed_Z_samps{k, i}) - con_lofi.State_Solve(oed_Z_samps{k, i});

        rand_Z_samps{k, i} = md_oed.Generate_Random_Design(N);
        rand_D_samps{k, i} = con_hifi.State_Solve(rand_Z_samps{k, i}) - con_lofi.State_Solve(rand_Z_samps{k, i});
        subrand_Z_samps{k, i} = md_oed.Generate_Random_Design_from_Subspace(N);
        subrand_D_samps{k, i} = con_hifi.State_Solve(subrand_Z_samps{k, i}) - con_lofi.State_Solve(subrand_Z_samps{k, i});
    end
end

save('OED_Ensemble_Results.mat');
