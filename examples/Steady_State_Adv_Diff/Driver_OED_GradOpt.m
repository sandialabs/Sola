% Clear Workspace and Add Interfaces to Path
OED_Setup;

% Perform Offline OED Computations
oed_interface = MD_OED_Interface_Diff(data_interface, con_lofi);
oed_reg_coeff = 10;
beta_0 = randn(num_evals, 1);

% Initialize Quantities
N = 5;
Jhat_seq_oed = zeros(N, 1);
oed_z_error = zeros(N, 1);
Z = [];
D = [];
betas = [];
z_bar = z_lofi;

% Sequential OED
md_oed = MD_OED_OptGrad(opt_prob_interface, data_interface, u_prior_interface, z_prior_interface, md_hessian_analysis, oed_interface);

for p = 1:N
    fprintf('\nStep %d:\n-------------\n', p);

    % Set Parameters for OED
    if p == 1
        z_p = z_lofi;
    else
        md_oed.Offline_Computation(md_cont_update, z_bar, u_bar);
        [beta_new, z_p] = md_oed.Generate_Seq_Optimal_Design(beta_0, alpha_d, oed_reg_coeff, betas, betas_cont(:, end));
        betas = [betas; beta_new];
        z_p = z_p(:, end);
        % disp(norm(z_p - z_bar)/norm(z_bar))
        % disp(norm(z_lofi - z_bar)/norm(z_bar))
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
    theta_post = Extract_mean_theta(md_post_sampling.post_data);

    % Obtain Optimal Solution Update via Continuation
    num_continuation_steps = 1;
    md_cont_update = MD_Continuation_Update(md_post_sampling, md_hessian_analysis, num_continuation_steps);
    [u_cont, z_cont, betas_cont] = md_cont_update.Posterior_Update_Mean_PC_beta();
    z_bar = z_cont(:, end);

    [delta_mean, ~] = md_post_sampling.Posterior_Discrepancy_Samples(z_bar);
    u_bar = con_lofi.State_Solve(z_bar) + delta_mean{:};

    % Obtain Optimal Solution Update via HDSA (linearization)
    % md_update = MD_Update(md_post_sampling, md_hessian_analysis);
    % z_bar_2 = md_update.Posterior_Update_Mean();
    % disp(norm(z_bar-z_bar_2)/norm(z_bar))
    % z_bar = z_lofi - PHinvB(theta_post);

    % Display Stats
    Jhat_seq_oed(p) = opt_hifi.Jhat(z_bar);
    oed_z_error(p) = oed_z_error_fn(z_bar);
    fprintf('Objective of z_bar: \t%.3f\n', Jhat_seq_oed(p));
    if p == 1
        fprintf('Percent Improvement: \t%.2f%%\n\n', 100 * (Jhat_lofi - Jhat_seq_oed(p)) / (Jhat_lofi - Jhat_hifi));
    else
        fprintf('Percent Improvement: \t%.2f%%\n\n', 100 * (Jhat_seq_oed(p - 1) - Jhat_seq_oed(p)) / (Jhat_seq_oed(p - 1) - Jhat_hifi));
    end

end

% Plot Objective Function over N
if true
    figure;
    hold on;
    xlim([0 N]);
    yline(Jhat_hifi, "k--", "DisplayName", "Hi-Fi", "LineWidth", 3, "Layer", "Bottom", "Alpha", 1);
    yline(Jhat_lofi, "r--", "DisplayName", "Lo-Fi", "LineWidth", 3, "Layer", "Bottom", "Alpha", 1);
    % yline(Jhat_best_proj, "b--", "DisplayName", "Best-Proj", "LineWidth", 3, "Layer", "Bottom", "Alpha", 1);
    plot(0:N, [Jhat_lofi; Jhat_seq_oed], ".-", "Color", "#00C83A", "DisplayName", "Sequential OED");
    xlabel("Evaluations ($N$)", "Interpreter", "latex");
    ylabel("Objective $\hat{J}(\cdot)$", "Interpreter", "latex");
    legend("location", "east", "Interpreter", "latex");
    title("Optimization Objective over Evals");
end
