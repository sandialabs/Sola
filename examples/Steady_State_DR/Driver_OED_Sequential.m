% Clear Workspace and Add Interfaces to Path
% clear;
close all;
% clc;
addpath(genpath('../../src'));
rng(0);

% Set Default Font Axes and Line Width
set(0, "DefaultAxesFontSize", 20);
set(0, "DefaultLineLineWidth", 3);
set(0, "DefaultLineMarkerSize", 20);

% Retrieve Model Parameters (D, Z, diff/reg/react_coeff, m, u_lofi, z_hifi/lofi; remove Z and D though)
load Optimization_Results.mat;
clear Z D;

% Set Hi-Fi and Lo-Fi Objectives and Constraints
obj = Diff_React_Objective(m, reg_coeff);
con_lofi = Diff_React_Constraint(m, diff_coeff, react_coeff);
opt_lofi = Reduced_Space_Optimization(obj, con_lofi);
con_hifi = Diff_React_HiFi_Constraint(con_lofi);
opt_hifi = Reduced_Space_Optimization(obj, con_hifi);
x = con_lofi.x;

% Show initial objective
fprintf("\nStep 0:\n-------------");
Jhat_lofi = opt_hifi.Jhat(z_lofi);
Jhat_hifi = opt_hifi.Jhat(z_hifi);
fprintf('Objective of z_lofi: \t%.3f\n', Jhat_lofi);
fprintf('Objective of z_hifi: \t%.3f\n\n', Jhat_hifi);

% Set Data Interface (no data there yet, except for z_lofi/u_lofi)
data_interface = MD_Data_Interface_Diff_React(u_lofi, z_lofi);

% Generate Priors for u and z
alpha_u = 2^2;
alpha_z = 1.e-10;
alpha_d = 1.e-4;
u_prior_interface = MD_Elliptic_u_Prior_Interface_Diff_React(alpha_u, opt_lofi);
z_prior_interface = MD_Elliptic_z_Prior_Interface_Diff_React(alpha_z, opt_lofi);

% Error with z_hifi
oed_z_error_fn = @(z) sqrt((z - z_hifi)' * z_prior_interface.Apply_M_z(z - z_hifi)) / sqrt(z_hifi' * z_prior_interface.Apply_M_z(z_hifi));

% %% Hessian analysis
opt_prob_interface = MD_Opt_Prob_Interface_Sabl(opt_lofi, data_interface);
md_hessian_analysis = MD_Hessian_Analysis(opt_prob_interface, z_prior_interface);
num_evals = 4;
oversampling = 20;
md_hessian_analysis.Compute_Hessian_GEVP(data_interface.z_opt, num_evals, oversampling);

% Perform Offline OED Computations - This is used to generate many random designs (to avoid OED in next steps)
alpha_zd = 1.e-2;
beta_zd = 1.e-2;
reg_coeff = 1.e-6;
beta_0 = randn(num_evals, 1);
data_interface = MD_Data_Interface_Diff_React(u_lofi, z_lofi);
oed_interface = MD_OED_Interface_Diff_React(con_lofi, alpha_zd, beta_zd);
% oed_interface = MD_OED_Interface_Diff_React(data_interface, con_lofi, alpha_zd, beta_zd);
% md_oed = MD_OED(opt_prob_interface, data_interface, u_prior_interface, z_prior_interface, md_hessian_analysis, oed_interface);
% md_oed.Offline_Computation();
% Z_random_samples = md_oed.Generate_Random_Design(100);

figure;
hold on;
plot(x, con_hifi.State_Solve(z_lofi), "r-", "DisplayName", "$S(\tilde{z})$");
plot(x, con_hifi.State_Solve(z_hifi), "k--", "DisplayName", "$S(z^*)$");
% plot(x, obj.T, "DisplayName", "Target");
title("Seq-OED State (Iteration 0)");
ylim([10 17]);
legend("Location", "northwest", "interpreter", "latex");
% saveas(gcf, "SeqOED_N_0.png")
%% Iterate for each data point
N = 6;
Jhat_oed = zeros(N, 1);
oed_z_error = zeros(N, 1);
z_bar = z_lofi;
Z = [];
D = [];
betas = [];
for p = 1:N
    % Update Data Interface (BEWARE: the change affects the dependencies directly!)
    fprintf('\nStep %d:\n-------------\n', p);

    % Sequential OED
    md_oed = MD_OED_Seq(opt_prob_interface, u_lofi, z_lofi, z_bar, u_prior_interface, z_prior_interface, md_hessian_analysis, oed_interface);
    md_oed.Offline_Computation();

    % Set Parameters for OED
    if p == 1
        z_p = z_lofi;
    else
        [beta_new, z_p] = md_oed.Generate_Seq_Optimal_Design(beta_0, alpha_d, reg_coeff, betas);
        betas = [betas; beta_new];
        z_p = z_p(:, end);
    end

    % Obtain Discrepancies
    % Z = Z_random_samples(:, 1:p);
    % Z = [Z z_bar+(Z_random_samples(:, p)-z_lofi)];
    % Z = [Z z_bar];
    Z = [Z z_p];
    D_p = Evaluate_Discrepancy(con_hifi, con_lofi, z_p);
    D = [D D_p];
    data_interface.Set_Z_and_D(Z, D);

    % Perform Posterior Sampling
    data_interface.z_opt = z_bar;
    md_post_sampling = MD_Posterior_Sampling(data_interface, u_prior_interface, z_prior_interface);
    md_post_sampling.Compute_Posterior_Data(alpha_d, 1);

    % Obtain Optimal Solution Update
    data_interface.z_opt = z_lofi;
    md_update = MD_Update(md_post_sampling, md_hessian_analysis);
    z_bar = md_update.Posterior_Update_Mean();

    % Display Stats
    Jhat_oed(p) = opt_hifi.Jhat(z_bar);
    oed_z_error(p) = oed_z_error_fn(z_bar);
    fprintf('Objective of z_bar: \t%.3f\n', Jhat_oed(p));
    % fprintf('Rel. Err of z_bar: \t%.2f%%\n', 100 * oed_z_error(p));
    % fprintf('Diff. w/ z_hifi obj.: \t%.2f%%\n', 100 * (Jhat_oed(p) - Jhat_hifi) / (Jhat_hifi));
    if p == 1
        fprintf('Percent Improvement: \t%.2f%%\n\n', 100 * (Jhat_lofi - Jhat_oed(p)) / (Jhat_lofi - Jhat_hifi));
    else
        fprintf('Percent Improvement: \t%.2f%%\n\n', 100 * (Jhat_oed(p - 1) - Jhat_oed(p)) / (Jhat_oed(p - 1) - Jhat_hifi));
    end

    % figure;
    % hold on;
    % plot(x, z_lofi, "r-", "DisplayName", "Lo-Fi Sol.");
    % plot(x, z_hifi, "k--", "DisplayName", "Hi-Fi Sol.");
    % plot(x, z_bar, "b-", "DisplayName", "Updated Sol.");
    % title("Lo-Fi & Hi-Fi Controls (Iteration " + p + ")");
    % legend("Location", "best");

    figure;
    hold on;
    plot(x, con_hifi.State_Solve(z_lofi), "r-", "DisplayName", "$S(\tilde{z})$");
    plot(x, con_hifi.State_Solve(z_hifi), "k--", "DisplayName", "$S(z^*)$");
    plot(x, con_hifi.State_Solve(z_bar), "b-", "DisplayName", "$S(\bar{z})$");
    ylim([10 17]);
    % plot(x, obj.T, "DisplayName", "Target");
    title("Seq-OED State (Iteration " + p + ")");
    legend("Location", "northwest", "interpreter", "latex");
    % saveas(gcf, "SeqOED_N_"+p+".png")

end

% figure;
% hold on;
% xlim([0 N])
% yline(Jhat_hifi, "k--", "DisplayName", "Hi-Fi", "LineWidth", 3, "Layer", "Bottom", "Alpha", 1)
% yline(Jhat_lofi, "r--", "DisplayName", "Lo-Fi", "LineWidth", 3, "Layer", "Bottom", "Alpha", 1)
% % plot(0:N, [Jhat_lofi; old_oed(1:N)], ".-", "Color", "#1F618D", "DisplayName", "Standard OED")
% plot(0:N, [Jhat_lofi; Jhat_oed], ".-", "Color", "#00C83A", "DisplayName", "Sequential OED")
% xlabel("Evaluations ($N$)", "Interpreter", "latex")
% ylabel("Objective $\hat{J}(\cdot)$", "Interpreter", "latex")
% legend("location", "east", "Interpreter", "latex")
% title("Optimization Objective over Evals")
