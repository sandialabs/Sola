%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Setup Synthetic Interfaces

%%
clear;
close all;
rng(121234);

m = 51;
n = m;
x = linspace(0, 1, m)';

data_interface = MD_Data_Interface_synthetic_test_OED(m);
data_interface.Load_Data();

u_prior_interface = MD_u_Prior_Interface_synthetic_test_OED(m);
z_prior_interface = MD_z_Prior_Interface_synthetic_test_OED(m);
md_prior_sampling = MD_Prior_Sampling(data_interface, u_prior_interface, z_prior_interface);

md_post_sampling = MD_Posterior_Sampling(data_interface, u_prior_interface, z_prior_interface);
alpha_d = 1.e-5;

opt_prob_interface = MD_Opt_Prob_Interface_synthetic_test_OED(m);
md_hessian_analysis = MD_Hessian_Analysis(opt_prob_interface, z_prior_interface);

num_evals = 30;
oversampling = 10;
md_hessian_analysis.Compute_Hessian_GEVP(data_interface.z_opt, num_evals, oversampling);
M_z_norm = @(z) sqrt(z' * z_prior_interface.Apply_M_z(z));
% Get the OED Setup ready
md_oed = MD_OED(opt_prob_interface, data_interface, u_prior_interface, z_prior_interface, md_hessian_analysis);
md_oed.Offline_Computation();

%% Perform OED (Reduced iterations for regression test)
N = 5;
Z = [];
D = [];
betas = [];
z_lofi = data_interface.z_opt;
Jhat_DC_oed = zeros(N, 1);
z_bars = zeros(n, N);
beta_bars = zeros(num_evals, N);
beta_0 = randn(num_evals, 1);
alpha_k_denom = trace(z_prior_interface.Apply_W_z_Inverse(z_prior_interface.M));

for p = 1:N
    % Set Parameters for OED
    if p == 1
        z_p = z_lofi;
    else
        if p == 2
            alpha_k_num = M_z_norm(z_bar - z_lofi)^2;
        else
            alpha_k_num = M_z_norm(z_bar - z_bars(:, p - 2))^2;
        end
        alpha_k = alpha_k_num / alpha_k_denom;
        md_oed.Set_Covariance_Coefficient(alpha_k);
        [beta_new, z_p] = md_oed.Generate_Seq_Optimal_Design(beta_0, alpha_d, betas, beta_bar, alpha_k_num);
        betas = [betas; beta_new];
        z_p = z_p(:, end);
    end

    % Obtain Discrepancies
    Z = [Z z_p];
    D_p = data_interface.Evaluate_Discrepancy(z_p);
    D = [D D_p];
    data_interface.Set_Z_and_D(Z, D);

    % Perform Posterior Sampling (TODO: Reuse data)
    md_post_sampling = MD_Posterior_Sampling(data_interface, u_prior_interface, z_prior_interface);
    md_post_sampling.Compute_Posterior_Data(alpha_d, 1);

    % Obtain Optimal Solution Update via Continuation
    num_continuation_steps = 3;
    md_cont_update = MD_Continuation_Update(md_post_sampling, md_hessian_analysis, num_continuation_steps);
    [u_cont, z_cont, betas_cont] = md_cont_update.Posterior_Update_Mean();
    z_bar = z_cont(:, end);
    z_bars(:, p) = z_bar;
    beta_bar = betas_cont(:, end);
    beta_bars(:, p) = beta_bar;

end

% --- Verification ---
ref_sol = load('reference_solution.mat');
rel_err = @(val_est, val_true) norm(val_est - val_true) / norm(val_true);
ref_diff = rel_err(Z, ref_sol.Z);

if ref_diff > 1.e-9
    fprintf(2, '\nmodel_discrepancy/synthetic_test_OED Test 1 failed.\n');
else
    fprintf(1, '\nmodel_discrepancy/synthetic_test_OED Test 1 passed.\n');
end
