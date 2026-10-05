%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clear Workspace and Add Interfaces to Path
addpath(genpath('..'));
addpath(genpath('../../../src'));
rng(0);

% Set Default Font Axes and Line Width
set(0, 'DefaultAxesFontSize', 20);
set(0, 'DefaultLineLineWidth', 3);
set(0, 'DefaultLineMarkerSize', 20);

% Retrieve model parameters (saved by Driver_Opt)
load Optimization_Results.mat;
clear Z D;
n = length(z_lofi);

% Set Hi-Fi and Lo-Fi objectives and constraints
obj = Subsurface_Objective(m, reg_coeff, p0);
con_lofi = Subsurface_LoFi_Constraint(m, k0, alpha, p0, viscosity);
opt_lofi = Reduced_Space_Optimization(obj, con_lofi);
con_hifi = Subsurface_HiFi_Constraint(con_lofi);
opt_hifi = Reduced_Space_Optimization(obj, con_hifi);
x = con_lofi.x;

% Compute objectives at z_hifi and z_lofi
Jhat_lofi = opt_hifi.Jhat(z_lofi);
Jhat_hifi = opt_hifi.Jhat(z_hifi);

% Note this does not contain access to Z/D until set explicitly.
data_interface = MD_Data_Interface_Subsurface(u_lofi, z_lofi);

% Generate priors for pressure discrepancy and control updates
alpha_u = 1;
alpha_z = 1.e-3;
alpha_d = (1.e-2)^2 * alpha_u;
u_prior_interface = MD_Elliptic_u_Prior_Interface_Subsurface(alpha_u, opt_lofi);
z_prior_interface = MD_Elliptic_z_Prior_Interface_Subsurface(alpha_z, opt_lofi);

% Convenience norms / diagnostics
M_z_norm = @(z) sqrt(z' * z_prior_interface.Apply_M_z(z));
W_z_norm = @(z) sqrt(z' * z_prior_interface.Apply_W_z(z));
oed_z_error_fn = @(z) M_z_norm(z - z_hifi) / M_z_norm(z_hifi);

% Hessian analysis
opt_prob_interface = MD_Opt_Prob_Interface_Sola(opt_lofi, data_interface);
md_hessian_analysis = MD_Hessian_Analysis(opt_prob_interface, z_prior_interface);
num_evals = 4;
oversampling = 20;
md_hessian_analysis.Compute_Hessian_GEVP(data_interface.z_opt, num_evals, oversampling);

% Get best possible z under projected problem
z_best_proj = z_lofi + md_hessian_analysis.evecs * (md_hessian_analysis.evecs \ (z_hifi - z_lofi));
Jhat_best_proj = opt_hifi.Jhat(z_best_proj);

% Display initial objectives
fprintf('\nStep 0:\n-------------\n');
fprintf('Objective of z_lofi: \t%.3f\n', 100 * Jhat_lofi);
fprintf('Objective of z_hifi: \t%.3f\n', 100 * Jhat_hifi);
fprintf('Objective of z_proj: \t%.3f\n\n', 100 * Jhat_best_proj);
