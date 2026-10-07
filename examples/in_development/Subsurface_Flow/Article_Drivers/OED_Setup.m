%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Clear Workspace and Add Interfaces to Path.  Build paths relative to
% this setup file so the driver can be launched from either this directory
% or the parent Subsurface_Flow directory.
article_dir = fileparts(mfilename('fullpath'));
example_dir = fileparts(article_dir);
repo_root = fullfile(example_dir, '..', '..', '..');
addpath(genpath(example_dir));
addpath(genpath(fullfile(repo_root, 'src')));
rng(0);

% Set Default Font Axes and Line Width
set(0, 'DefaultAxesFontSize', 20);
set(0, 'DefaultLineLineWidth', 3);
set(0, 'DefaultLineMarkerSize', 20);

% Retrieve model parameters (saved by Driver_Opt).  Prefer a local copy
% for article workflows, but fall back to the parent example directory.
results_file = fullfile(article_dir, 'Optimization_Results.mat');
if ~isfile(results_file)
    results_file = fullfile(example_dir, 'Optimization_Results.mat');
end
load(results_file);
if ~exist('Z', 'var') || ~exist('D', 'var')
    error(['Optimization_Results.mat must contain discrepancy design data Z and D. ', ...
           'Run examples/in_development/Subsurface_Flow/Driver_Opt.m again.']);
end
hyperparam_Z = Z;
hyperparam_D = D;

% Guard against silently using optimization data from the earlier weak-
% discrepancy version of the example.
expected_alpha = 0.5;
expected_reg_coeff = 1.e-6;
expected_hifi_leakoff_coeff = 30;
expected_hifi_leakoff_center = 0.55;
expected_hifi_leakoff_width = 0.25;
expected_lofi_permeability_model = 'exponential';
if ~exist('lofi_permeability_model', 'var') || ...
        ~strcmp(lofi_permeability_model, expected_lofi_permeability_model) || ...
        ~exist('hifi_leakoff_coeff', 'var') || ...
        ~exist('hifi_leakoff_center', 'var') || ...
        ~exist('hifi_leakoff_width', 'var') || ...
        abs(alpha - expected_alpha) > 1.e-12 || ...
        abs(reg_coeff - expected_reg_coeff) > 1.e-15 || ...
        abs(hifi_leakoff_coeff - expected_hifi_leakoff_coeff) > 1.e-12 || ...
        abs(hifi_leakoff_center - expected_hifi_leakoff_center) > 1.e-12 || ...
        abs(hifi_leakoff_width - expected_hifi_leakoff_width) > 1.e-12
    error(['Optimization_Results.mat was generated with old subsurface parameters. ', ...
           'Run examples/in_development/Subsurface_Flow/Driver_Opt.m again before the article drivers.']);
end

n = length(z_lofi);

% Set Hi-Fi and Lo-Fi objectives and constraints
obj = Subsurface_Objective(m, reg_coeff, p0);
con_lofi = Subsurface_LoFi_Constraint(m, k0, alpha, p0, viscosity);
opt_lofi = Reduced_Space_Optimization(obj, con_lofi);
con_hifi = Subsurface_HiFi_Constraint(con_lofi, hifi_leakoff_coeff, hifi_leakoff_center, hifi_leakoff_width);
opt_hifi = Reduced_Space_Optimization(obj, con_hifi);
x = con_lofi.x;

% Compute objectives at z_hifi and z_lofi
Jhat_lofi = opt_hifi.Jhat(z_lofi);
Jhat_hifi = opt_hifi.Jhat(z_hifi);

% Data interface used for automated prior hyperparameter initialization.
% Driver_Opt saves a small set of discrepancy evaluations; these are used here to
% estimate correlation lengths, prior magnitudes, and the data-noise scale.
data_interface = MD_Data_Interface_Subsurface(u_lofi, z_lofi);
data_interface.Set_Z_and_D(hyperparam_Z, hyperparam_D);

% Automatically determine discrepancy-prior hyperparameters from the saved
% discrepancy data.  This replaces the old manually specified alpha_u,
% alpha_z, and alpha_d values while using the common numeric Laplacian prior
% infrastructure in src/model_discrepancy.
u_hyperparam_interface = MD_u_Hyperparameter_Interface_Subsurface(x, false);
u_prior_interface = MD_Numeric_Laplacian_u_Prior_Interface(con_lofi.S, con_lofi.M, data_interface, u_hyperparam_interface);
u_hyperparam_interface.alpha_d = 1e3*u_hyperparam_interface.alpha_d;


z_hyperparam_interface = MD_z_Hyperparameter_Interface_Subsurface(x, con_lofi);
z_prior_interface = MD_Numeric_Laplacian_z_Prior_Interface(con_lofi.S, con_lofi.M, data_interface, z_hyperparam_interface, u_prior_interface);

alpha_u = u_prior_interface.alpha_u;
alpha_z = z_prior_interface.alpha_z;
alpha_d = u_hyperparam_interface.alpha_d;
beta_u = u_prior_interface.beta_u;
beta_z = z_prior_interface.beta_z;

fprintf('Automated discrepancy hyperparameters:\n');
fprintf('  alpha_u = %.4e, beta_u = %.4e, alpha_d = %.4e\n', alpha_u, beta_u, alpha_d);
fprintf('  alpha_z = %.4e, beta_z = %.4e\n', alpha_z, beta_z);

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
