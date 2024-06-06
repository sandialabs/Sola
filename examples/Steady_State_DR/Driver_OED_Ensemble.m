% Clear Workspace and Add Interfaces to Path
clear;
close all;
% clc;
addpath(genpath('../../src'));

% Set Default Font Axes and Line Width
set(0, "DefaultAxesFontSize", 20);
set(0, "DefaultLineLineWidth", 3);
set(0, "DefaultLineMarkerSize", 20);

% Retrieve Model Parameters (D, Z, diff/reg/react_coeff, m, u_lofi, z_hifi/lofi)
load Optimization_Results.mat;

% Set Hi-Fi and Lo-Fi Objectives and Constraints
obj = Diff_React_Objective(m, reg_coeff);
con_lofi = Diff_React_Constraint(m, diff_coeff, react_coeff);
opt_lofi = Reduced_Space_Optimization(obj, con_lofi);
con_hifi = Diff_React_HiFi_Constraint(con_lofi);
opt_hifi = Reduced_Space_Optimization(obj, con_hifi);
x = con_lofi.x;

% Load Data (NOTE: FROM Optimization_Results.mat) - PROCEED WITH CAUTION ON ANYTHING USING DATA_INTERFACE!
data_interface = MD_Data_Interface_Diff_React();
data_interface.Load_Data();

% Generate Priors for u and z
alpha_u = 2^2;
alpha_z = 1.e-10;
u_prior_interface = MD_Elliptic_u_Prior_Interface_Diff_React(alpha_u, opt_lofi);
z_prior_interface = MD_Elliptic_z_Prior_Interface_Diff_React(alpha_z, opt_lofi);

% %% Hessian analysis - Only uses z_opt from data_interface
opt_prob_interface = MD_Opt_Prob_Interface_Sabl(opt_lofi, data_interface);
md_hessian_analysis = MD_Hessian_Analysis(opt_prob_interface, z_prior_interface);
num_evals = 4;
oversampling = 20;
md_hessian_analysis.Compute_Hessian_GEVP(data_interface.z_opt, num_evals, oversampling);

% Perform Offline OED Computations - USES data_interface
alpha_zd = 0.05;
beta_zd = 1.e-2;
oed_interface = MD_OED_Interface_Diff_React(data_interface, con_lofi, alpha_zd, beta_zd);
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

% k dictates what N is
for k = 1:p
    N = N_range(k);
    for i = 1:samps_per_N
        beta_0 = randn(num_evals * (N - 1), 1);
        alpha_d = 1.e-4;
        reg_coeff = 1.e-6;
        [oed_beta_samps{k, i}, oed_Z_samps{k, i}] = md_oed.Generate_Optimal_Design(beta_0, alpha_d, reg_coeff);
        oed_obj(k, i) = md_oed.Evaluate_OED_Objective(oed_beta_samps{k, i}, alpha_d, reg_coeff);
        rand_Z_samps{k, i} = md_oed.Generate_Random_Design(N);
        subrand_Z_samps{k, i} = md_oed.Generate_Random_Design_from_Subspace(N);

        % Note that the PDE solver can only handle one state at a time.
        for ix = 1:size(oed_Z_samps{k, i}, 2)
            oed_D_samps{k, i}(:, ix) = con_hifi.State_Solve(oed_Z_samps{k, i}(:, ix)) - con_lofi.State_Solve(oed_Z_samps{k, i}(:, ix));
            rand_D_samps{k, i}(:, ix) = con_hifi.State_Solve(rand_Z_samps{k, i}(:, ix)) - con_lofi.State_Solve(rand_Z_samps{k, i}(:, ix));
            subrand_D_samps{k, i}(:, ix) = con_hifi.State_Solve(subrand_Z_samps{k, i}(:, ix)) - con_lofi.State_Solve(subrand_Z_samps{k, i}(:, ix));
        end

    end
end

save('OED_Ensemble_Results.mat');
disp("Saved.");
