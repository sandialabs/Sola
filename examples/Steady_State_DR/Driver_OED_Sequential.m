% Get the OED Setup ready
OED_Setup;
oed_reg_coeff = 1.e-3;
beta_0 = randn(num_evals, 1);

% Iterate over steps
N = 5;
Jhat_seq_oed = zeros(N, 1);
oed_z_error = zeros(N, 1);
Z = [];
D = [];
betas = [];
z_bar = z_lofi;

% Saving quantities
% seq_oed_mean_theta = cell(N, 1);
% seq_oed_mean_z = cell(N, 1);
% seq_oed_Z = cell(N, 1);

% Sequential OED
md_oed = MD_OED_Seq(opt_prob_interface, data_interface, u_prior_interface, z_prior_interface, md_hessian_analysis, oed_interface);
md_oed.Offline_Computation();

for p = 1:N
    % Update Data Interface (with prior center)
    fprintf('\nStep %d:\n-------------\n', p);
    data_interface.Update_z_opt(z_bar);

    % Set Parameters for OED
    if p == 1
        z_p = z_lofi;
    else
        % beta_0 = beta_cont(:, end);
        % if p == 2
        %     oed_reg_coeffs = [1.e-6, 1.e-5, 1.e-4, 2.e-4, 5.e-4, 1.e-3, 1.e-2];
        %     [beta_lc, Z_lc, post_var, reg_val] = md_oed.L_Curve_Analysis(beta_0, alpha_d, oed_reg_coeffs, betas);
        % end
        [beta_new, z_p] = md_oed.Generate_Seq_Optimal_Design(beta_0, alpha_d, oed_reg_coeff, betas);
        betas = [betas; beta_new];
        z_p = z_p(:, end);
        % disp(norm(z_p - z_bar))
        % z_p = z_bar;
    end

    % Obtain Discrepancies
    Z = [Z z_p];
    D_p = Evaluate_Discrepancy(con_hifi, con_lofi, z_p);
    D = [D D_p];
    data_interface.Set_Z_and_D(Z, D);
    seq_oed_Z{p} = Z;

    % Perform Posterior Sampling (TODO: Reuse data)
    md_post_sampling = MD_Posterior_Sampling(data_interface, u_prior_interface, z_prior_interface);
    md_post_sampling.Compute_Posterior_Data(alpha_d, 1);

    % Obtain Optimal Solution Update via Continuation
    num_continuation_steps = 2;
    md_cont_update = MD_Continuation_Update(md_post_sampling, md_hessian_analysis, num_continuation_steps);
    [u_cont, z_cont, beta_cont] = md_cont_update.Posterior_Update_Mean_PC_beta();
    z_bar = z_cont(:, end);

    % Obtain Optimal Solution Update via HDSA (linearization)
    % md_update = MD_Update(md_post_sampling, md_hessian_analysis);
    % z_bar_2 = md_update.Posterior_Update_Mean();
    % disp(norm(z_bar-z_bar_2)/norm(z_bar))

    % Display Stats
    Jhat_seq_oed(p) = opt_hifi.Jhat(z_bar);
    oed_z_error(p) = oed_z_error_fn(z_bar);
    fprintf('Objective of z_bar: \t%.3f\n', Jhat_seq_oed(p));
    if p == 1
        fprintf('Percent Improvement: \t%.2f%%\n\n', 100 * (Jhat_lofi - Jhat_seq_oed(p)) / (Jhat_lofi - Jhat_hifi));
    else
        fprintf('Percent Improvement: \t%.2f%%\n\n', 100 * (Jhat_seq_oed(p - 1) - Jhat_seq_oed(p)) / (Jhat_seq_oed(p - 1) - Jhat_hifi));
    end

    % seq_oed_mean_theta{p} = Extract_mean_theta(md_post_sampling.post_data);
    % seq_oed_mean_z{p} = z_bar;

end

% Plot Objective Function over N
if true
    figure;
    hold on;
    xlim([0 N]);
    yline(Jhat_hifi, "k--", "DisplayName", "Hi-Fi", "LineWidth", 3, "Layer", "Bottom", "Alpha", 1);
    yline(Jhat_lofi, "r--", "DisplayName", "Lo-Fi", "LineWidth", 3, "Layer", "Bottom", "Alpha", 1);
    plot(0:N, [Jhat_lofi; Jhat_seq_oed], ".-", "Color", "#00C83A", "DisplayName", "Sequential OED");
    xlabel("Evaluations ($N$)", "Interpreter", "latex");
    ylabel("Objective $\hat{J}(\cdot)$", "Interpreter", "latex");
    legend("location", "east", "Interpreter", "latex");
    title("Optimization Objective over Evals");
end

% save('Seq_OED_Results.mat', 'Jhat_seq_oed', 'seq_oed_mean_theta', 'seq_oed_mean_z', 'seq_oed_Z');
