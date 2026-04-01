% Get the OED Setup ready
OED_Setup;

% Perform Offline OED Computations - USES data_interface
md_oed = MD_OED(opt_prob_interface, data_interface, u_prior_interface, z_prior_interface, md_hessian_analysis);
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
alpha_k_denom = trace(z_prior_interface.Apply_W_z_Inverse(z_prior_interface.M));

for p = 1:N
    fprintf('\nStep %d:\n-------------\n', p);

    % Set Parameters for OED
    if p == 1
        z_p = z_lofi;
    else
        tdisp("Starting OED");
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
    tdisp("Evaluating discrepancy");
    Z = [Z z_p];
    D_p = Evaluate_Discrepancy(con_hifi, con_lofi, max(z_p, 0));
    D = [D D_p];
    data_interface.Set_Z_and_D(Z, D);
    tdisp("Discrepancy evaluated");

    % Perform Posterior Sampling (TODO: Reuse data)
    tdisp("Start posterior computation");
    md_post_sampling = MD_Posterior_Sampling(data_interface, u_prior_interface, z_prior_interface);
    md_post_sampling.Compute_Posterior_Data(alpha_d, 1);
    tdisp("End posterior computation");

    % Obtain Optimal Solution Update via Continuation
    num_continuation_steps = 3;
    md_cont_update = MD_Continuation_Update(md_post_sampling, md_hessian_analysis, num_continuation_steps);
    tdisp("Start continuation");
    [u_cont, z_cont, betas_cont] = md_cont_update.Posterior_Update_Mean();
    tdisp("End continuation");
    z_bar = z_cont(:, end);
    z_bars(:, p) = z_bar;
    beta_bar = betas_cont(:, end);
    beta_bars(:, p) = beta_bar;

    % Display Stats
    Jhat_DC_oed(p) = Jhat_hifi_fn(z_bar);
    fprintf('Objective of z_bar: \t%.3f\n', Jhat_DC_oed(p));
    if p == 1
        fprintf('Percent Improvement: \t%.2f%%\n\n', 100 * (Jhat_lofi - Jhat_DC_oed(p)) / (Jhat_lofi - Jhat_hifi));
    else
        fprintf('Percent Improvement: \t%.2f%%\n\n', 100 * (Jhat_DC_oed(p - 1) - Jhat_DC_oed(p)) / (Jhat_DC_oed(p - 1) - Jhat_hifi));
    end

end

% Plot Objective Function over N;
if true
    figure;
    hold on;
    xlim([0 5]);
    yline(1e-4 * Jhat_hifi, "k--", "DisplayName", "Hi-Fi", "LineWidth", 3, "Layer", "Bottom", "Alpha", 1);
    yline(1e-4 * Jhat_lofi, "r--", "DisplayName", "Lo-Fi", "LineWidth", 3, "Layer", "Bottom", "Alpha", 1);
    plot(0:5, 1e-4 * [Jhat_lofi; Jhat_DC_oed(1:5)], ".-", "Color", "#BAB86C", "DisplayName", "Solution Updates");
    xlabel("Evaluations ($N$)", "Interpreter", "latex");
    ylabel("High-fidelity objective", "Interpreter", "latex");
    legend("location", "east", "Interpreter", "latex");
    saveas(gca, 'NS-1', 'epsc');
    % title("Optimization Objective over Evals");
end

% figure;
% semilogy(0:5, 1e-4 * ([Jhat_lofi; Jhat_DC_oed(1:5)] - Jhat_hifi), ".-", "Color", "#BAB86C", "DisplayName", "Solution Updates");
% hold on;
% xlim([0 5]);
% % yline(1e-4 * Jhat_hifi, "k--", "DisplayName", "Hi-Fi", "LineWidth", 3, "Layer", "Bottom", "Alpha", 1);
% yline(1e-4 * Jhat_lofi, "r--", "DisplayName", "Lo-Fi", "LineWidth", 3, "Layer", "Bottom", "Alpha", 1);
% xlabel("Evaluations ($N$)", "Interpreter", "latex");
% ylabel("Objective $\hat{J}(\cdot)$", "Interpreter", "latex");
% legend("location", "east", "Interpreter", "latex");

Z_oed = Z;
D_oed = D;

figure;
hold on;
plot(flip(z_lofi), "Color", [0.7 0.7 0.7], "DisplayName", "$\tilde{z}$", "LineWidth", 3);
plot(flip(z_hifi), "k-", "DisplayName", "$z_{\rm hifi}$", "LineWidth", 3);
plot(flip(z_bars(:, 5)), "b-", "DisplayName", "$\overline{z}_5$", "LineWidth", 3);
xlabel("$x$", "Interpreter", "latex");
ylabel("$c(x, 0)$", "Interpreter", "latex");
legend("location", "northeast", "Interpreter", "latex");
saveas(gca, 'NS-2', 'epsc');
