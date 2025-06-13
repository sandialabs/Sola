%% Set up
addpath(genpath('../../src'));
load Optimization_Results.mat;

obj = Adv_Diff_Objective(m, reg_coeff);
con_hifi = Adv_Diff_Constraint(m, diff_coeff, vel_coeff, robin_coeff);
con_lofi = Diff_Constraint(con_hifi);
opt_hifi = Reduced_Space_Optimization(obj, con_hifi);
opt_lofi = Reduced_Space_Optimization(obj, con_lofi);
x = con_lofi.x;

warning("Before running this script, ensure that hyperparameters (e.g., alpha_*, oed_reg_coeff, num_evals) match.");
%% HDSA interfaces
data_interface = MD_Data_Interface_Diff();
data_interface.Load_Data();

alpha_u = 15;
alpha_z = 5;
alpha_d = (1.e-2)^2 * alpha_u;
u_prior_interface = MD_Elliptic_u_Prior_Interface_Diff(alpha_u, opt_lofi);
z_prior_interface = MD_Elliptic_z_Prior_Interface_Diff(alpha_z, opt_lofi);
oed_reg_coeff = 1.e-6;

%% Hessian analysis
opt_prob_interface = MD_Opt_Prob_Interface_Sabl(opt_lofi, data_interface);
md_hessian_analysis = MD_Hessian_Analysis(opt_prob_interface, z_prior_interface);
num_evals = 4;
oversampling = 20;
md_hessian_analysis.Compute_Hessian_GEVP(data_interface.z_opt, num_evals, oversampling);

%% OED
oed_interface = MD_OED_Interface_Diff(data_interface, con_lofi);
md_oed = MD_OED(opt_prob_interface, data_interface, u_prior_interface, z_prior_interface, md_hessian_analysis, oed_interface);
md_oed.verbosity = false;
md_oed.Offline_Computation();
N = 6;

oed_beta_samps = cell(N, 1);
oed_Z_samps = cell(N, 1);
oed_D_samps = cell(N, 1);
oed_obj = zeros(N, 1);

oed_Z_samps{1, 1} = z_lofi;
oed_beta_samps{1, 1} = nan;
oed_obj(1, 1) = nan;
for k = 1:N
    if k ~= 1
        beta_0 = randn(num_evals * (k - 1), 1);
        [oed_beta_samps{k, 1}, oed_Z_samps{k, 1}] = md_oed.Generate_Optimal_Design(beta_0, alpha_d, oed_reg_coeff);
        oed_obj(k, 1) = md_oed.Evaluate_OED_Objective(oed_beta_samps{k, 1}, alpha_d, oed_reg_coeff);
    end
    oed_D_samps{k, 1} = con_hifi.State_Solve(oed_Z_samps{k, 1}) - con_lofi.State_Solve(oed_Z_samps{k, 1});
end

save('OED_Ensemble_Results.mat');
