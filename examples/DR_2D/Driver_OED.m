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

alpha_u = 2^2;
alpha_z = 1.e-4; % 1.e-6;
u_prior_interface = MD_Elliptic_u_Prior_Interface_Diff_React(alpha_u, opt);
z_prior_interface = MD_Elliptic_z_Prior_Interface_Diff_React(alpha_z, opt);

%% Hessian analysis
opt_prob_interface = MD_Opt_Prob_Interface_Sabl(opt, data_interface);
md_hessian_analysis = MD_Hessian_Analysis(opt_prob_interface, z_prior_interface);
num_evals = 8; % 4;
oversampling = 20;
md_hessian_analysis.Compute_Hessian_GEVP(data_interface.z_opt, num_evals, oversampling);

%% OED
oed_interface = MD_OED_Interface_Diff_React(data_interface, con);

md_oed = MD_OED(opt_prob_interface, data_interface, u_prior_interface, z_prior_interface, md_hessian_analysis, oed_interface);
md_oed.Offline_Computation();

N = 4;
% Z = md_oed.Generate_Random_Design(N);
% for k = 1:N
%    mesh.Plot_Field(Z(:,k),['Random sample ',num2str(k)]);
%    set(gca, 'fontsize', 18);
% end
%
% Z = md_oed.Generate_Random_Design_from_Subspace(N);
% for k = 1:N
%    mesh.Plot_Field(Z(:,k),['Random subspace sample ',num2str(k)]);
%    set(gca, 'fontsize', 18);
% end

beta_0 = 10 * randn(num_evals * (N - 1), 1);
alpha_d = 1.e-2;
% reg_coeffs = 10.^(-6:-1:-10)';
% [beta_L_curve, Z_L_curve, post_var, reg_val] = md_oed.L_Curve_Analysis(beta_0, alpha_d, reg_coeffs);
% figure;
% plot(post_var, reg_val, 'o', 'MarkerSize', 10);
% set(gca, 'fontsize', 18);
% L_curve_index = 3;
% Z = Z_L_curve{L_curve_index};

reg_coeff = 1.e-8;
[beta, Z] = md_oed.Generate_Optimal_Design(beta_0, alpha_d, reg_coeff);

mesh.Plot_Field(Z(:, 1), 'OED sample 1');
set(gca, 'fontsize', 18);
for k = 2:N
    mesh.Plot_Field(Z(:, k) - Z(:, 1), ['OED diff sample ', num2str(k)]);
    set(gca, 'fontsize', 18);
end
