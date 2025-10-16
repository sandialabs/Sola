% Clear Workspace and Add Interfaces to Path
addpath(genpath('../../src'));
rng(0);

% Set Default Font Axes and Line Width
set(0, "DefaultAxesFontSize", 20);
set(0, "DefaultLineLineWidth", 3);
set(0, "DefaultLineMarkerSize", 20);

% Set Python environment and variables
pyenv('Version', '/usr/local/anaconda3/envs/FenicsEnvCompat/bin/python', 'ExecutionMode', 'InProcess');
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
alpha_z = 1.e-1; % If too large, may cause bouncing.
alpha_u = (1)^2; % smallar alpha_u -> smaller expected magnitude of discrepancy
alpha_d = (1.e-2)^2 * alpha_u; % smaller alpha_d -> more certainty (linearity) in data
beta_t = 50;
z_prior_interface = MD_Elliptic_z_Prior_Interface_Tracer(alpha_z, opt_lofi);

% Set Transient Prior
n_t = 25;
n_y = 31;
T = 0.1;
u_hyperparam_interface = MD_u_Hyperparameter_Interface(true);
u_hyperparam_interface.Set_beta_t(beta_t);
u_hyperparam_interface.Set_alpha_u(alpha_u);
u_hyperparam_interface.Set_alpha_d(alpha_d);
spatial_u_prior_interface = MD_Elliptic_u_Prior_Interface_Tracer(alpha_u, opt_lofi);
transient_prior_cov = MD_Transient_Prior_Covariance_Sabl(data_interface, u_hyperparam_interface, T, n_t, n_y);
u_prior_interface = MD_Transient_Elliptic_u_Prior_Interface(data_interface, spatial_u_prior_interface, transient_prior_cov);

% Calculate Relative OED Error (Lambda Function for now)
M_z_norm = @(z) sqrt(z' * z_prior_interface.Apply_M_z(z));
W_z_norm = @(z) sqrt(z' * z_prior_interface.Apply_W_z(z));
oed_z_error_fn = @(z) M_z_norm(z - z_hifi) / M_z_norm(z_hifi);

% Hessian analysis
opt_prob_interface = MD_Opt_Prob_Interface_Python(data_interface, opt_lofi);
md_hessian_analysis = MD_Hessian_Analysis(opt_prob_interface, z_prior_interface);
num_evals = 4;
oversampling = 5;
disp("Computing Hessian GEVP...");
md_hessian_analysis.Compute_Hessian_GEVP(data_interface.z_init, num_evals, oversampling);

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
