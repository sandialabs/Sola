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

% Retrieve Model Parameters (D, Z, diff/reg/react_coeff, m, u_lofi, z_hifi/lofi)
load Optimization_Results.mat;
n = m;

% Set Hi-Fi and Lo-Fi Objectives and Constraints
obj = Adv_Diff_Objective(m, reg_coeff);
con_hifi = Adv_Diff_Constraint(m, diff_coeff, vel_coeff, robin_coeff);
con_lofi = Diff_Constraint(con_hifi);
opt_hifi = Reduced_Space_Optimization(obj, con_hifi);
opt_lofi = Reduced_Space_Optimization(obj, con_lofi);
x = con_lofi.x;

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
md_oed = MD_OED(opt_prob_interface, data_interface, u_prior_interface, z_prior_interface, md_hessian_analysis, oed_interface);
md_oed.Offline_Computation();

% Set Parameters for OED
N = 5;
rng(0);
beta_0 = randn(num_evals * (N - 1), 1);
reg_coeff = 1.e-6;
[betas, Z] = md_oed.Generate_Optimal_Design(beta_0, alpha_d, reg_coeff);

% Generate Design (Generate_Random_Design(N), Generate_Random_Design_from_Subspace(N), Generate_Optimal_Design(...))
% Z = md_oed.Generate_Random_Design(N);
D = Evaluate_Discrepancy(con_hifi, con_lofi, Z);
data_interface.Set_Z_and_D(Z, D);

% % Sample from Posterior (i.e., solve problem) - USES DATA INTERFACE
num_post_samples = 1;
md_post_sampling = MD_Posterior_Sampling(data_interface, u_prior_interface, z_prior_interface);
md_post_sampling.Compute_Posterior_Data(alpha_d, num_post_samples);
% [delta_mean, delta_samples] = md_post_sampling.Posterior_Discrepancy_Samples(Z);

% % Sample from Posterior of Optimal Solution
md_update = MD_Update(md_post_sampling, md_hessian_analysis);
[z_update_mean, z_update_samples] = md_update.Posterior_Update_Samples();

% Sample from Design Prior of Z
num_prior_samples = 1;
Z_prior_samps = md_oed.Generate_Random_Design(num_prior_samples);

% Evaluate Objectives
fprintf("\n\n");
fprintf("\nObjective Value of Hi-Fi Control: \t" + opt_hifi.Jhat(z_hifi));
fprintf("\nObjective Value of Lo-Fi Control: \t" + opt_hifi.Jhat(z_lofi));
fprintf("\nObjective Value of Updated Control: \t" + opt_hifi.Jhat(z_update_mean));
fprintf("\n\n");
fprintf("\nError of Lo-Fi Control: \t" + oed_z_error_fn(z_lofi));
fprintf("\nError of Updated Control: \t" + oed_z_error_fn(z_update_mean));

if false
    % Comparison of Lo-Fi and Hi-Fi Control Solutions
    figure;
    hold on;
    plot(x, Z_prior_samps, "Color", [0.9, 0.9, 0.7], "HandleVisibility", "off");
    plot(x, z_update_samples, "Color", [0.7, 0.9, 0.9], "HandleVisibility", "off");
    plot(x, z_lofi, "Color", 0.5 * [0.9, 0.9, 0.3], "DisplayName", "Lo-Fi Sol.");
    plot(x, z_hifi, "DisplayName", "Hi-Fi Sol.");
    plot(x, z_update_mean, "Color", 0.5 * [0.3, 0.9, 0.9], "DisplayName", "Updated Sol.");
    title("Lo-Fi & Hi-Fi Controls");
    legend("Location", "best");
    % z_prior_samples_w = z_prior_interface.Sample_with_sCovariance_W_z_Inverse(10);
    % figure;
    % plot(x, z_prior_samples)

    % Comparison of Lo-Fi and Hi-Fi State Solutions
    figure;
    hold on;
    plot(x, con_hifi.State_Solve(z_lofi), "r-", "DisplayName", "Lo-Fi Sol.");
    plot(x, con_hifi.State_Solve(z_hifi), "k--", "DisplayName", "Hi-Fi Sol.");
    plot(x, con_hifi.State_Solve(z_update_mean), "b-", "DisplayName", "Updated Sol.");
    % plot(x, obj.T, "DisplayName", "Target");
    title("States for Lo-Fi & Hi-Fi Controls");
    legend("Location", "best");

    figure;
    hold on;
    plot(x, z_lofi, "r-", "DisplayName", "Lo-Fi Sol.");
    plot(x, z_hifi, "k--", "DisplayName", "Hi-Fi Sol.");
    plot(x, z_update_mean, "b-", "DisplayName", "Updated Sol.");
    % plot(x, obj.T, "DisplayName", "Target");
    title("Lo-Fi & Hi-Fi Controls");
    legend("Location", "best");
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kronecker product computation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Im = eye(m);
Mz = z_prior_interface.M;
B1 = @(x) opt_prob_interface.Apply_Solution_Operator_z_Jacobian_Transpose(opt_prob_interface.Apply_Misfit_Hessian([Im kron(Im, z_lofi' * Mz)] * x, u_lofi, z_lofi), z_lofi);
B2 = @(x) [zeros(n, n) kron(opt_prob_interface.Misfit_Gradient(u_lofi, z_lofi)', Mz)] * x;
B = @(x) B1(x) + B2(x);
PHinvB = @(x) md_hessian_analysis.Apply_Projected_RS_Hessian_Inverse(B(x));

%% Making sure PHinvB works by applying it on the posterior mean...
% theta_post_mean_tmp = md_update.Posterior_Theta_Mean_Temp();
% z_update_tmp = z_lofi - PHinvB(theta_post_mean_tmp);
% fprintf("\n\n")
% disp(norm(z_update_mean - z_update_tmp)/norm(z_update_mean))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Obtaining Best Theta
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
A_lofi = con_hifi.diff_coeff * con_hifi.S + con_hifi.robin_coeff * con_hifi.robin_bc;
discrep_mat = 10^2 * (inv(A_lofi + con_hifi.vel_coeff * con_hifi.V) - inv(A_lofi));
best_theta = [zeros(m, 1); reshape(discrep_mat', m * n, 1)];

%% Making sure best_theta is correct by evaluating the discrepancy at different points
% eval_discrep_theta = @(z, theta) [Im kron(Im, z' * Mz)] * theta;
% discrep_zero = con_hifi.State_Solve(z_hifi)-con_lofi.State_Solve(z_hifi);
% discrep_est = eval_discrep_theta(z_hifi, best_theta);
% disp(norm(discrep_zero-discrep_est)/norm(discrep_zero))

z_best_HDSA = z_lofi - PHinvB(best_theta);
Jhat_HDSA = opt_hifi.Jhat(z_best_HDSA);
fprintf("\n");
fprintf("\nError of Best-HDSA Control: \t" + oed_z_error_fn(z_best_HDSA)); % 0.14764
fprintf("\nObjective of Best-HDSA Control: \t" + Jhat_HDSA); % 19.9269

best_z = z_best_HDSA;
save('Truth_Results.mat', 'best_theta', 'best_z', 'con_hifi', 'con_lofi', 'opt_hifi', 'opt_lofi', 'obj', 'z_hifi', 'z_lofi', 'u_prior_interface', 'z_prior_interface', 'opt_prob_interface', 'md_hessian_analysis');
