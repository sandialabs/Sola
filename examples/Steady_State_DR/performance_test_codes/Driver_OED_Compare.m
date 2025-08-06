% Get the OED Setup ready
addpath(genpath('..'));
addpath(genpath('../../../src'));
OED_Setup;

% Perform Offline OED Computations - USES data_interface
md_oed = MD_OED_DeltaCov(opt_prob_interface, data_interface, u_prior_interface, z_prior_interface, md_hessian_analysis, oed_interface);
md_oed.Offline_Computation();

%% Perform OED
N = 5;
Z = [];
D = [];
betas = [];
Jhat_DC_oed = zeros(N, 1);
oed_reg_coeff = 1.e-2;
z_bars = zeros(n, N);
beta_0 = randn(num_evals, 1);

%% Step 1
p = 1;
fprintf('\nStep %d:\n-------------\n', p);
z_p = z_lofi;
Z = [Z z_p];
D_p = Evaluate_Discrepancy(con_hifi, con_lofi, z_p);
D = [D D_p];
data_interface.Set_Z_and_D(Z, D);

% Perform Posterior Sampling
md_post_sampling = MD_Posterior_Sampling(data_interface, u_prior_interface, z_prior_interface);
md_post_sampling.Compute_Posterior_Data(alpha_d, 1);

% Obtain Optimal Solution Update
num_continuation_steps = 3;
md_cont_update = MD_Continuation_Update(md_post_sampling, md_hessian_analysis, num_continuation_steps);
[u_cont, z_cont, betas_cont] = md_cont_update.Posterior_Update_Mean_PC_beta();
z_bar = z_cont(:, end);
z_bars(:, p) = z_bar;
fprintf('Objective at z_bar: \t%.3f\n', opt_hifi.Jhat(z_bar));

%% Step 2 - OED
p = 2;
fprintf('\nStep %d (OED):\n-------------\n', p);

covar_coeff = W_z_norm(z_bar - z_lofi)^2 / n;
% disp(sqrt(covar_coeff));
md_oed.Set_Covariance_Coefficient(covar_coeff);
% [beta_l, Z_l, post_var_l, reg_val_l] = md_oed.L_Curve_Analysis(beta_0, alpha_d, [0, 10.^(-4:0.5:-1)], betas, betas_cont(:, end));
% md_oed.Evaluate_Posterior_Cov_Trace(beta_0, alpha_d, betas_cont(:, end));
% beta_0 = betas_cont(:, end);

beta_0 = randn(num_evals * 2, 1);
[beta_new, z_p] = md_oed.Generate_Seq_Optimal_Design(beta_0, alpha_d, oed_reg_coeff, betas, betas_cont(:, end));
disp(size(z_p));
betas = [betas; beta_new];
disp(norm(z_p - z_bar) / norm(z_lofi - z_bar));

% Obtain Discrepancies
Z = [Z z_p];
D_p = Evaluate_Discrepancy(con_hifi, con_lofi, z_p);
D = [D D_p];
data_interface.Set_Z_and_D(Z, D);

% Perform Posterior Sampling
md_post_sampling = MD_Posterior_Sampling(data_interface, u_prior_interface, z_prior_interface);
md_post_sampling.Compute_Posterior_Data(alpha_d, 1);

% Obtain Optimal Solution Update
num_continuation_steps = 3;
md_cont_update = MD_Continuation_Update(md_post_sampling, md_hessian_analysis, num_continuation_steps);
[u_cont, z_cont, betas_cont] = md_cont_update.Posterior_Update_Mean_PC_beta();
z_bar_OED = z_cont(:, end);
z_bars(:, p) = z_bar_OED;
fprintf('Objective at z_bar: \t%.3f\n', opt_hifi.Jhat(z_bar_OED));

% %% Step 2 - Random
% p = 2;
% num_samples = 2;
% objectives = zeros(num_samples, 1);
% covar_coeff = W_z_norm(z_bar - z_lofi)^2 / n;
% % z_ps = z_bar + sqrt(covar_coeff) * z_prior_interface.Sample_with_Covariance_W_z_Inverse(num_samples);
% z_ps = z_bar + sqrt(covar_coeff) * md_hessian_analysis.evecs * (md_hessian_analysis.evecs' * z_prior_interface.Sample_with_Covariance_W_z_Inverse(num_samples))/W_z_norm(z_lofi);
% fprintf('\nStep %d (Rand.):\n-------------\n', p);
% rng(2);

% wb = waitbar(0, 'Starting');

% for i = 1:num_samples
%     waitbar(i / num_samples, wb, sprintf('Progress: %d %%', floor(i / num_samples * 100)));
%     z_p = z_ps(:, i);

%     % Obtain Discrepancies
%     Z = [Z(:, 1:p - 1) z_p];
%     D_p = Evaluate_Discrepancy(con_hifi, con_lofi, z_p);
%     D = [D(:, 1:p - 1) D_p];
%     data_interface.Set_Z_and_D(Z, D);

%     % Perform Posterior Sampling
%     md_post_sampling = MD_Posterior_Sampling(data_interface, u_prior_interface, z_prior_interface);
%     md_post_sampling.Compute_Posterior_Data(alpha_d, 1);

%     % Obtain Optimal Solution Update
%     num_continuation_steps = 3;
%     md_cont_update = MD_Continuation_Update(md_post_sampling, md_hessian_analysis, num_continuation_steps);
%     [u_cont, z_cont, betas_cont] = md_cont_update.Posterior_Update_Mean_PC_beta();
%     z_bar = z_cont(:, end);

%     % Compute and store the objective
%     objectives(i) = opt_hifi.Jhat(z_bar);
%     fprintf('Objective at z_bar (Sample %d): \t%.3f\n', i, objectives(i));
% end

% % Display all objectives
% fprintf('\nAll Objectives:\n');
% disp(objectives);

% covar_coeff = W_z_norm(z_bar - z_bars(:, p - 2))^2 / n;
%% Plot Objective Function over N
% if true
% figure;
% hold on;
% xlim([0 N]);
% yline(Jhat_hifi, "k--", "DisplayName", "Hi-Fi", "LineWidth", 3, "Layer", "Bottom", "Alpha", 1);
% yline(Jhat_lofi, "r--", "DisplayName", "Lo-Fi", "LineWidth", 3, "Layer", "Bottom", "Alpha", 1);
% plot(0:N, [Jhat_lofi; Jhat_DC_oed], ".-", "Color", "#BAB86C", "DisplayName", "DeltaCov OED");
% xlabel("Evaluations ($N$)", "Interpreter", "latex");
% ylabel("Objective $\hat{J}(\cdot)$", "Interpreter", "latex");
% legend("location", "east", "Interpreter", "latex");
% title("Optimization Objective over Evals");
% end
