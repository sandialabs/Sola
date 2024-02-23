%% Set up
clear;
close all;
clc;
addpath(genpath('../../src'));
load Optimization_Results.mat;

mesh = PDE_Meshing(h);
x = mesh.x;
y = mesh.y;
M = mesh.M;
diff_react_lofi = Diff_React_Lofi(mesh);
diff_react_hifi = Diff_React_Hifi(mesh);
m = length(x);

reg_coeff = 1.e-4;
obj = Diff_React_Objective(diff_react_lofi, reg_coeff);
con = Diff_React_Constraint(diff_react_lofi);
opt = Reduced_Space_Optimization(obj, con);

%% HDSA interfaces
data_interface = MD_Data_Interface_Diff_React();
data_interface.Load_Data();
opt_prob_interface = MD_Opt_Prob_Interface_Sabl(opt, data_interface);
alpha_u = 2^2;
alpha_z = (1 / 50000)^2;
u_prior_interface = MD_Elliptic_u_Prior_Interface_Diff_React(alpha_u, opt);
z_prior_interface = MD_Elliptic_z_Prior_Interface_Diff_React(alpha_z, opt);

%% Hessian analysis
md_hessian_analysis = MD_Hessian_Analysis(opt_prob_interface, z_prior_interface);
num_evals = 4;
oversampling = 20;
md_hessian_analysis.Compute_Hessian_GEVP(data_interface.z_opt, num_evals, oversampling);

%% OED
oed_interface = MD_OED_Interface_Diff_React(data_interface, con);

md_oed = MD_OED(opt_prob_interface, data_interface, u_prior_interface, z_prior_interface, md_hessian_analysis, oed_interface);
md_oed.Offline_Computation();

samps_per_N = 5;
N_range = (2:4)';
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
        beta_0 = 10 * randn(num_evals * (N - 1), 1);
        alpha_d = 1.e-2;
        reg_coeff = 1.e-8;

        [oed_beta_samps{k, i}, oed_Z_samps{k, i}] = md_oed.Generate_Optimal_Design(beta_0, alpha_d, reg_coeff);
        rand_Z_samps{k, i} = md_oed.Generate_Random_Design(N);
        subrand_Z_samps{k, i} = md_oed.Generate_Random_Design_from_Subspace(N);

        oed_D_samps{k, i} = zeros(m, N);
        rand_D_samps{k, i} = zeros(m, N);
        subrand_D_samps{k, i} = zeros(m, N);
        for j = 1:N
            tmp = diff_react_hifi.State_Solve(diff_react_hifi.Map_z_to_Control_Fun(oed_Z_samps{k, i}(:, j)));
            oed_D_samps{k, i}(:, j) = tmp - con.State_Solve(oed_Z_samps{k, i}(:, j));

            tmp = diff_react_hifi.State_Solve(diff_react_hifi.Map_z_to_Control_Fun(rand_Z_samps{k, i}(:, j)));
            rand_D_samps{k, i}(:, j) = tmp - con.State_Solve(rand_Z_samps{k, i}(:, j));

            % tmp = diff_react_hifi.State_Solve(diff_react_hifi.Map_z_to_Control_Fun(subrand_Z_samps{k,i}(:,j)));
            % subrand_D_samps{k, i}(:,j) = tmp - con.State_Solve(subrand_Z_samps{k, i}(:,j));
        end
    end
end

save('OED_Ensemble_Results.mat');
