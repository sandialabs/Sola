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
u_lofi = con_lofi.State_Solve(z_lofi);
% u_lofi = load("data/lofi_optim_sol.mat").k_opt_lofi;
z_hifi = load("data/hifi_optim_sol.mat").k0_hifi;

% Show initial objective
fprintf("\nStep 0:\n-------------");
Jhat_lofi = Jhat_hifi_fn(z_lofi);
Jhat_hifi = Jhat_hifi_fn(z_hifi);
fprintf('Objective of z_lofi: \t%.3f\n', Jhat_lofi);
fprintf('Objective of z_hifi: \t%.3f\n\n', Jhat_hifi);

% Set Data Interface
data_interface = MD_Data_Interface_Tracer(u_lofi, z_lofi);

% Generate Priors for u and z
alpha_d = 1.e-7; % Controls speed of improvement (large)
alpha_z = 0.1; % If too large, may cause bouncing.
alpha_u = (0.2)^2; % Larger prior if large.
beta_t = 50;
beta_i = 1.e5; % Minimal Impact
num_evals = 6; % significant impact when very small; very little beyond
oversampling = 1;
z_prior_interface = MD_Elliptic_z_Prior_Interface_Tracer(alpha_z, opt_lofi);

% Set Transient Prior
n_t = 25;
n_y = 31;
T = 0.1;
transient_prior_cov = MD_Transient_Prior_Covariance_Sabl(beta_t, beta_i, T, n_t, n_y);
u_prior_interface = MD_Transient_Elliptic_u_Prior_Interface_Tracer(alpha_u, transient_prior_cov, opt_lofi);

% num_prior_samples = 100;
% md_prior_sampling = MD_Prior_Sampling(data_interface, u_prior_interface, z_prior_interface);
% delta_samples = md_prior_sampling.Prior_Discrepancy_Samples_at_z_opt(num_prior_samples);
% plot(x, delta_samples(end-30:end, :))

% Error with z_hifi
oed_z_error_fn = @(z) sqrt((z - z_hifi)' * z_prior_interface.Apply_M_z(z - z_hifi)) / sqrt(z_hifi' * z_prior_interface.Apply_M_z(z_hifi));

% Perform Hessian Analysis
opt_prob_interface = MD_Opt_Prob_Interface_Python(data_interface);
md_hessian_analysis = MD_Hessian_Analysis(opt_prob_interface, z_prior_interface);
disp("Computing Hessian GEVP...");

md_hessian_analysis.Compute_Hessian_GEVP(data_interface.z_init, num_evals, oversampling);

% Perform Offline OED Computations
alpha_zd = 1.e-1;
beta_zd = 1.e-1;
reg_coeff = 1.e-12;
beta_0 = randn(num_evals, 1);
oed_interface = MD_OED_Interface_Tracer(data_interface, con_lofi, alpha_zd, beta_zd);

% Plot low-fidelity and high-fidelity states
% pyplot(x, con_hifi.State_Solve(z_lofi), 'r-', x, con_hifi.State_Solve(z_hifi), 'k--', 'Legend', {'Low-Fidelity', 'High-Fidelity'});

%% Iterate for each data point
N = 5;
Jhat_oed = zeros(N, 1);
oed_z_error = zeros(N, 1);
Z = [];
D = [];
betas = [];
z_bar = z_lofi;

disp("Discrep. Eval. for Hifi...");
D_z_hifi = Evaluate_Discrepancy(con_hifi, con_lofi, z_hifi);
disp("Discrep. Eval. for Hifi [DONE].");

for p = 1:N
    % Update Data Interface (with prior center)
    fprintf('\nStep %d:\n-------------\n', p);
    data_interface.Update_z_opt(z_bar);

    % Sequential OED
    md_oed = MD_OED_Seq(opt_prob_interface, data_interface, u_prior_interface, z_prior_interface, md_hessian_analysis, oed_interface);
    md_oed.Offline_Computation();

    % Set Parameters for OED
    if p == 1
        z_p = z_lofi;
        % betas = [betas; 0*beta_0]; % This is for the updated sequential OED
    else
        % z_p = z_bar;
        [beta_new, z_p] = md_oed.Generate_Seq_Optimal_Design(beta_0, alpha_d, reg_coeff, betas);
        betas = [betas; beta_new];
        z_p = z_p(:, end); % Redundancy for standard OED.
    end

    % Obtain Discrepancies
    Z = [Z z_p];
    D_p = Evaluate_Discrepancy(con_hifi, con_lofi, z_p);
    D = [D D_p];

    % Update data_interface
    data_interface.Set_Z_and_D(Z, D);
    data_interface.Update_z_opt(z_bar);

    % Perform Posterior Sampling (TODO: Avoid recomputing computed data points)
    md_post_sampling = MD_Posterior_Sampling(data_interface, u_prior_interface, z_prior_interface);
    md_post_sampling.Compute_Posterior_Data(alpha_d, 1);
    % post_mean_z_hifi(:, p) = cell2mat(md_post_sampling.Posterior_Discrepancy_Samples(z_hifi));
    % post_mean_z_p(:, p) = cell2mat(md_post_sampling.Posterior_Discrepancy_Samples(z_p));
    % disp("Discrepancy Norm at Hifi z");
    % disp(vecnorm(post_mean_z_hifi(:, p)));
    % disp("Discrepancy Norm at Data z");
    % disp(vecnorm(post_mean_z_p(:, p)));

    % Obtain Optimal Solution Update
    md_update = MD_Update(md_post_sampling, md_hessian_analysis);
    z_bar = md_update.Posterior_Update_Mean();

    % Display Stats
    Jhat_oed(p) = Jhat_hifi_fn(z_bar);
    fprintf('Objective of z_bar: \t%.3f\n', Jhat_oed(p));
    if p == 1
        fprintf('Percent Improvement: \t%.2f%%\n\n', 100 * (Jhat_lofi - Jhat_oed(p)) / (Jhat_lofi - Jhat_hifi));
    else
        fprintf('Percent Improvement: \t%.2f%%\n\n', 100 * (Jhat_oed(p - 1) - Jhat_oed(p)) / (Jhat_oed(p - 1) - Jhat_hifi));
    end
end

% pyplot(0:N, [Jhat_lofi; Jhat_oed], '.-', 'Title', 'Optimization Objective over Evals');
show_figures = true;
if show_figures
    figure;
    hold on;
    xlim([0 N]);
    yline(Jhat_hifi, "k--", "DisplayName", "Hi-Fi", "LineWidth", 3, "Layer", "Bottom", "Alpha", 1);
    yline(Jhat_lofi, "r--", "DisplayName", "Lo-Fi", "LineWidth", 3, "Layer", "Bottom", "Alpha", 1);
    % plot(0:N, [Jhat_lofi; old_oed(1:N)], ".-", "Color", "#1F618D", "DisplayName", "Standard OED")
    plot(0:N, [Jhat_lofi; Jhat_oed], ".-", "Color", "#00C83A", "DisplayName", "Sequential OED");
    xlabel("Evaluations ($N$)", "Interpreter", "latex");
    ylabel("Objective $\hat{J}(\cdot)$", "Interpreter", "latex");
    legend("location", "east", "Interpreter", "latex");
    title("Optimization Objective over Evals");
end
