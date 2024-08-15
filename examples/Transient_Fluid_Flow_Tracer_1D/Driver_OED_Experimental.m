% Clear Workspace and Add Interfaces to Path
addpath(genpath('../../src'));
rng(0);

% Set Default Font Axes and Line Width
set(0, "DefaultAxesFontSize", 20);
set(0, "DefaultLineLineWidth", 3);
set(0, "DefaultLineMarkerSize", 20);


% Set Python environment and variables
pyenv('Version', '/usr/local/anaconda3/envs/FenicsEnvCompat/bin/python'); 
setenv('PKG_CONFIG_PATH', '/usr/local/anaconda3/envs/FenicsEnvCompat/lib/pkgconfig');
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
% x = con_lofi.x;
disp("Breakpoint 1")

% Obtain high-fidelity and low-fidelity optimizersß
z_lofi = load("data/lofi_optim_sol.mat").k0_opt_lofi;
u_lofi = load("data/lofi_optim_sol.mat").k_opt_lofi;
z_hifi = load("data/hifi_optim_sol.mat").k0_hifi;


% Set Data Interface
data_interface = MD_Data_Interface_Tracer(u_lofi, z_lofi);

% Generate Priors for u and z
alpha_u = 2^2;
alpha_z = 1.e-10;
alpha_d = 1.e-4;
u_prior_interface = MD_Elliptic_u_Prior_Interface_Tracer(alpha_u, opt_lofi);
z_prior_interface = MD_Elliptic_z_Prior_Interface_Tracer(alpha_z, opt_lofi);
disp("Breakpoint 2")

% Error with z_hifi
oed_z_error_fn = @(z) sqrt((z - z_hifi)' * z_prior_interface.Apply_M_z(z - z_hifi)) / sqrt(z_hifi' * z_prior_interface.Apply_M_z(z_hifi));


% Perform Hessian Analysis
opt_prob_interface = MD_Opt_Prob_Interface_Sabl(opt_lofi, data_interface);
md_hessian_analysis = MD_Hessian_Analysis(opt_prob_interface, z_prior_interface);
num_evals = 4;
oversampling = 1;
disp("Breakpoint 3")

md_hessian_analysis.Compute_Hessian_GEVP(data_interface.z_init, num_evals, oversampling);



% TODO: Give access to x-values & compute vertex values for plotting


% 