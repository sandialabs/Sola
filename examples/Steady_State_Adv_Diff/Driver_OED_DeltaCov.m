% Clear Workspace and Add Interfaces to Path
% clear;
% close all;
% clc;
addpath(genpath('../../src'));
rng(0);

% Set Default Font Axes and Line Width
set(0, "DefaultAxesFontSize", 20);
set(0, "DefaultLineLineWidth", 3);
set(0, "DefaultLineMarkerSize", 20);

% Retrieve Model Parameters (D, Z, diff/reg/react_coeff, m, u_lofi, z_hifi/lofi)
load Optimization_Results.mat;

% Set Hi-Fi and Lo-Fi Objectives and Constraints
obj = Adv_Diff_Objective(m, reg_coeff);
con_hifi = Adv_Diff_Constraint(m, diff_coeff, vel_coeff, robin_coeff);
con_lofi = Diff_Constraint(con_hifi);
opt_hifi = Reduced_Space_Optimization(obj, con_hifi);
opt_lofi = Reduced_Space_Optimization(obj, con_lofi);
x = con_lofi.x;

% Show initial objective
fprintf("\nStep 0:\n-------------\n");
Jhat_lofi = opt_hifi.Jhat(z_lofi);
Jhat_hifi = opt_hifi.Jhat(z_hifi);
fprintf('Objective of z_lofi: \t%.3f\n', Jhat_lofi);
fprintf('Objective of z_hifi: \t%.3f\n\n', Jhat_hifi);

% Note this doesn't contain access to Z/D yet.
data_interface = MD_Data_Interface_Diff(u_lofi, z_lofi);

% Generate Priors for u and z
alpha_u = 15;
alpha_z = 5;
alpha_d = (1.e-2)^2 * alpha_u;
u_prior_interface = MD_Elliptic_u_Prior_Interface_Diff(alpha_u, opt_lofi);
z_prior_interface = MD_Elliptic_z_Prior_Interface_Diff(alpha_z, opt_lofi);

% Calculate Relative OED Error (Lambda Function for now)
oed_z_error_fn = @(z) sqrt((z - z_hifi)' * z_prior_interface.Apply_M_z(z - z_hifi)) / sqrt(z_hifi' * z_prior_interface.Apply_M_z(z_hifi));

% %% Hessian analysis
opt_prob_interface = MD_Opt_Prob_Interface_Sabl(opt_lofi, data_interface);
md_hessian_analysis = MD_Hessian_Analysis(opt_prob_interface, z_prior_interface);
num_evals = 4;
oversampling = 20;
md_hessian_analysis.Compute_Hessian_GEVP(data_interface.z_opt, num_evals, oversampling);

% Perform Offline OED Computations - USES data_interface
oed_interface = MD_OED_Interface_Diff(data_interface, con_lofi);
md_oed = MD_OED_DeltaCov(opt_prob_interface, data_interface, u_prior_interface, z_prior_interface, md_hessian_analysis, oed_interface);
md_oed.Offline_Computation();

%% Perform OED
N = 5;
Z = [];
D = [];
betas = [];
Jhat_DC_oed = zeros(N, 1);
oed_reg_coeff = 1.e-3;
beta_0 = randn(num_evals, 1);

for p = 1:N
    % Update Data Interface (with prior center)
    fprintf('\nStep %d:\n-------------\n', p);

    % Set Parameters for OED
    if p == 1
        z_p = z_lofi;
    else
        [beta_new, z_p] = md_oed.Generate_Seq_Optimal_Design(beta_0, alpha_d, oed_reg_coeff, betas, betas_cont(:, end));
        betas = [betas; beta_new];
        z_p = z_p(:, end);
    end

    % Obtain Discrepancies
    Z = [Z z_p];
    D_p = Evaluate_Discrepancy(con_hifi, con_lofi, z_p);
    D = [D D_p];
    data_interface.Set_Z_and_D(Z, D);

    % Perform Posterior Sampling (TODO: Reuse data)
    md_post_sampling = MD_Posterior_Sampling(data_interface, u_prior_interface, z_prior_interface);
    md_post_sampling.Compute_Posterior_Data(alpha_d, 1);
    theta_post = Extract_mean_theta(md_post_sampling.post_data);

    % Obtain Optimal Solution Update via Continuation
    num_continuation_steps = 1;
    md_cont_update = MD_Continuation_Update(md_post_sampling, md_hessian_analysis, num_continuation_steps);
    [u_cont, z_cont, betas_cont] = md_cont_update.Posterior_Update_Mean_PC_beta();
    z_bar = z_cont(:, end);

    % Display Stats
    Jhat_DC_oed(p) = opt_hifi.Jhat(z_bar);
    fprintf('Objective of z_bar: \t%.3f\n', Jhat_DC_oed(p));
    if p == 1
        fprintf('Percent Improvement: \t%.2f%%\n\n', 100 * (Jhat_lofi - Jhat_DC_oed(p)) / (Jhat_lofi - Jhat_hifi));
    else
        fprintf('Percent Improvement: \t%.2f%%\n\n', 100 * (Jhat_DC_oed(p - 1) - Jhat_DC_oed(p)) / (Jhat_DC_oed(p - 1) - Jhat_hifi));
    end

end

% Plot Objective Function over N
show_figures = true;
if show_figures
    figure;
    hold on;
    xlim([0 N]);
    yline(Jhat_hifi, "k--", "DisplayName", "Hi-Fi", "LineWidth", 3, "Layer", "Bottom", "Alpha", 1);
    yline(Jhat_lofi, "r--", "DisplayName", "Lo-Fi", "LineWidth", 3, "Layer", "Bottom", "Alpha", 1);
    plot(0:N, [Jhat_lofi; Jhat_DC_oed], ".-", "Color", "#BAB86C", "DisplayName", "DeltaCov OED");
    xlabel("Evaluations ($N$)", "Interpreter", "latex");
    ylabel("Objective $\hat{J}(\cdot)$", "Interpreter", "latex");
    legend("location", "east", "Interpreter", "latex");
    title("Optimization Objective over Evals");
end
