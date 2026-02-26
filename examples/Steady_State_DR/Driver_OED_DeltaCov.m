% Get the OED Setup ready
OED_Setup;

% Perform Offline OED Computations - USES data_interface
md_oed = MD_OED_DeltaCov(opt_prob_interface, data_interface, u_prior_interface, z_prior_interface, md_hessian_analysis);
md_oed.Offline_Computation();

%% Perform OED
N = 5;
Z = [];
D = [];
betas = [];
Jhat_DC_oed = zeros(N, 1);
z_bars = zeros(n, N);
beta_bars = zeros(num_evals, N);
beta_0 = randn(num_evals, 1);
% opt_hifi.Jhat(z_lofi + V*beta_bar)
alpha_k_denom = trace(z_prior_interface.Apply_W_z_Inverse(z_prior_interface.M));
alpha_k_denom_proj = trace(md_hessian_analysis.evecs' * z_prior_interface.M * md_hessian_analysis.evecs) / W_z_norm(z_lofi)^2;

for p = 1:N
    fprintf('\nStep %d:\n-------------\n', p);

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
    D_p = Evaluate_Discrepancy(con_hifi, con_lofi, z_p);
    D = [D D_p];
    data_interface.Set_Z_and_D(Z, D);

    % Perform Posterior Sampling (TODO: Reuse data)
    md_post_sampling = MD_Posterior_Sampling(data_interface, u_prior_interface, z_prior_interface);
    md_post_sampling.Compute_Posterior_Data(alpha_d, 1);

    % Obtain Optimal Solution Update via Continuation
    num_continuation_steps = 3;
    md_cont_update = MD_Continuation_Update(md_post_sampling, md_hessian_analysis, num_continuation_steps);
    [u_cont, z_cont, betas_cont] = md_cont_update.Posterior_Update_Mean_PC_beta();
    z_bar = z_cont(:, end);
    z_bars(:, p) = z_bar;
    beta_bar = betas_cont(:, end);
    beta_bars(:, p) = beta_bar;

    % Display Stats
    Jhat_DC_oed(p) = opt_hifi.Jhat(z_bar);
    fprintf('Objective of z_bar: \t%.5f\n', Jhat_DC_oed(p));
    if p == 1
        fprintf('Percent Improvement: \t%.2f%%\n\n', 100 * (Jhat_lofi - Jhat_DC_oed(p)) / (Jhat_lofi - Jhat_hifi));
    else
        fprintf('Percent Improvement: \t%.2f%%\n\n', 100 * (Jhat_DC_oed(p - 1) - Jhat_DC_oed(p)) / (Jhat_DC_oed(p - 1) - Jhat_hifi));
    end

end

% Plot Objective Function over N
if true
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

Z_oed = Z;
D_oed = D;
% save("performance_test_codes/oed-results.mat", "z_bars", "Jhat_DC_oed", "Z_oed", "D_oed");
