% Clear Workspace and Add Interfaces to Path
addpath(genpath('../../src'));
rng(0);

% Set Default Font Axes and Line Width
set(0, "DefaultAxesFontSize", 20);
set(0, "DefaultLineLineWidth", 3);
set(0, "DefaultLineMarkerSize", 20);

% Set Python environment and variables
pyenv('Version', '/usr/local/anaconda3/envs/FenicsEnvNew/bin/python', 'ExecutionMode', 'InProcess');
pythonFilePath = 'python';
if count(py.sys.path, pythonFilePath) == 0
    insert(py.sys.path, int32(0), pythonFilePath);
end

% Import and Reload the Python module
fluid_flow_1d_lofi = py.importlib.import_module("fluid_flow_1d_lofi");
py.importlib.reload(fluid_flow_1d_lofi);
fluid_flow_1d_hifi_eval = py.importlib.import_module("fluid_flow_1d_hifi_eval");
py.importlib.reload(fluid_flow_1d_hifi_eval);

% Set Hi-Fi and Lo-Fi Objectives and Constraints
obj = Tracer_Objective();
con_lofi = Tracer_LoFi_Constraint();
opt_lofi = Reduced_Space_Optimization(obj, con_lofi);
con_hifi = Tracer_HiFi_Constraint();
opt_hifi = Reduced_Space_Optimization(obj, con_hifi);
x = con_lofi.x;

% HIFI JHAT
Jhat_hifi_fn = @(z) obj.J(con_hifi.State_Solve_Terminal(z), z);
Jhat_lofi_fn = @(z) obj.J(con_lofi.State_Solve_Terminal(z), z);

% Obtain high-fidelity and low-fidelity optimizersß
z_lofi = load("data/lofi_optim_sol.mat").k0_opt_lofi;
n = length(z_lofi);
u_lofi = con_lofi.State_Solve(z_lofi);
disp("Evaluating Hi-Fi PDE at Hi-Fi Sol...");
z_hifi = load("data/hifi_optim_sol.mat").k0_hifi;
u_hifi = con_hifi.State_Solve(z_hifi);

% Show initial objective
fprintf("\nStep 0:\n-------------\n");
Jhat_lofi = Jhat_hifi_fn(z_lofi);
Jhat_hifi = Jhat_hifi_fn(z_hifi);
fprintf('Objective of z_lofi: \t%.3f\n', Jhat_lofi);
fprintf('Objective of z_hifi: \t%.3f\n\n', Jhat_hifi);

% Set Data Interface
data_interface = MD_Data_Interface_Tracer(u_lofi, z_lofi);

% Generate Priors for u and z
% alpha_u = 0*(1.e-5)^2;
% alpha_d = 0*(1.e-2)^2 * alpha_u;
% beta_t = 0*50;
alpha_z = 51.5;
alpha_u = 1.5277e-04;
alpha_d = 3.6663e-09;
beta_t = 0.0228;
% z_prior_interface = MD_Elliptic_z_Prior_Interface_Tracer(alpha_z, opt_lofi);
% determine_z_hyperparams = MD_Determine_z_Hyperparameters(data_interface, z_hyperparam_interface, u_prior_interface);
% determine_z_hyperparams.Determine_beta_z(); %0.0031
% z_prior_interface = MD_Elliptic_z_Prior_Interface_Tracer(1, opt_lofi);
% determine_z_hyperparams.Determine_alpha_z(z_prior_interface);
% Set Transient Prior
n_t = 25;
n_y = 31;
T = 0.1;
% u_hyperparam_interface = MD_u_Hyperparameter_Interface(true);
u_hyperparam_interface = MD_u_Hyperparameter_Interface_Tracer(n_t, n_y, true);

% u_prior_interface_new = MD_Numeric_Laplacian_u_Prior_Interface(pde_meshing.S, pde_meshing.M, data_interface, u_hyperparam_interface);

% u_hyperparam_interface.Set_beta_t(beta_t);
% u_hyperparam_interface.Set_alpha_u(alpha_u);
% u_hyperparam_interface.Set_alpha_d(alpha_d);
Z = load("results_great.mat", "Z_oed").Z_oed;
D = load("results_great.mat", "D_oed").D_oed;
data_interface.Set_Z_and_D(Z, D);

spatial_u_prior_interface = MD_Elliptic_u_Prior_Interface_Tracer(alpha_u, opt_lofi);
transient_prior_cov = MD_Transient_Prior_Covariance_Sabl(data_interface, u_hyperparam_interface, T, n_t, n_y);
u_prior_interface = MD_Transient_Elliptic_u_Prior_Interface(data_interface, spatial_u_prior_interface, transient_prior_cov);

z_hyperparam_interface = MD_z_Hyperparameter_Interface_Tracer(x);
z_prior_interface = MD_Numeric_Laplacian_z_Prior_Interface(con_lofi.S_z, con_lofi.M_z, data_interface, z_hyperparam_interface, u_prior_interface);

% Calculate Relative OED Error (Lambda Function for now)
M_z_norm = @(z) sqrt(z' * z_prior_interface.Apply_M_z(z));
W_z_norm = @(z) sqrt(z' * z_prior_interface.Apply_W_z(z));
oed_z_error_fn = @(z) M_z_norm(z - z_hifi) / M_z_norm(z_hifi);

% md_hessian_analysis.evecs' * z_prior_interface.Apply_W_z(Z_oed - z_lofi)/7.6251
tdisp = @(msg) fprintf('%s (%s)\n', msg, datetime('now'));

% Hessian analysis
opt_prob_interface = MD_Opt_Prob_Interface_Python(data_interface, opt_lofi);
md_hessian_analysis = MD_Hessian_Analysis(opt_prob_interface, z_prior_interface);
num_evals = 4;
oversampling = 5;
tdisp("Starting Hessian GEVP computation");
md_hessian_analysis.Compute_Hessian_GEVP(data_interface.z_init, num_evals, oversampling);
tdisp("Completed Hessian GEVP computation");

% Get OED Interface ready
alpha_zd = 1.e-1;
beta_zd = 1.e-1;
oed_interface = MD_OED_Interface_Tracer(data_interface, con_lofi, alpha_zd, beta_zd);

% Get best possible z under projected problem
z_best_proj = z_lofi + md_hessian_analysis.evecs * (md_hessian_analysis.evecs \ (z_hifi - z_lofi));
Jhat_best_proj = Jhat_hifi_fn(z_best_proj);

% Display initial objectives
fprintf("\nStep 0:\n-------------\n");
fprintf('Objective of z_lofi: \t%.3f\n', Jhat_lofi);
fprintf('Objective of z_hifi: \t%.3f\n', Jhat_hifi);
fprintf('Objective of z_proj: \t%.3f\n\n', Jhat_best_proj);

% Computing Hessian GEVP...

% Step 0:
% -------------
% Objective of z_lofi:    31.774
% Objective of z_hifi:    1.158
% Objective of z_proj:    8.750

% Step 1:
% -------------
% Objective of z_bar:     12.445
% Percent Improvement:    63.14%

% Step 2:
% -------------
% Objective of z_bar:     10.555
% Percent Improvement:    16.74%

% Step 3:
% -------------
% Objective of z_bar:     10.293
% Percent Improvement:    2.79%

% Step 4:
% -------------
% Objective of z_bar:     10.069
% Percent Improvement:    2.45%

% Step 5:
% -------------
% Objective of z_bar:     9.879
% Percent Improvement:    2.13%

% Best possible value w/ zero: 9.8724

% figure;
% hold on;
% shift_val = 9.8724
% xlim([0 N]);
% % yline(1.158, "k--", "DisplayName", "Hi-Fi", "LineWidth", 3, "Layer", "Bottom", "Alpha", 1);
% yline(31.774-shift_val, "r--", "DisplayName", "Lo-Fi", "LineWidth", 3, "Layer", "Bottom", "Alpha", 1);
% plot(0:N, [Jhat_lofi; 12.445; 10.555; 10.293; 10.069; 9.879]-shift_val, ".-", "Color", "#BAB86C", "DisplayName", "DeltaCov OED");
% xlabel("Evaluations ($N$)", "Interpreter", "latex");
% ylabel("Objective $\hat{J}(\cdot)$", "Interpreter", "latex");
% legend("location", "east", "Interpreter", "latex");
% title("Optimization Objective over Evals");
