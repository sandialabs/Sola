% Get the OED Setup ready
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
oed_reg_coeff = 0;
z_bars = zeros(n, N);
beta_bars = zeros(num_evals, N);
beta_0 = randn(num_evals, 1);

for p = 1:N
    fprintf('\nStep %d:\n-------------\n', p);

    % Set Parameters for OED
    if p == 1
        z_p = z_lofi;
    elseif p == 2
        z_p = z_hifi;
    else
        if p == 2
            covar_coeff = W_z_norm(z_bar - z_lofi)^2 / n;
            delta_beta = abs(beta_bar);
        else
            covar_coeff = W_z_norm(z_bar - z_bars(:, p - 2))^2 / n;
            delta_beta = abs(beta_bar - beta_bars(:, p - 2));
        end

        nonlcon = @(beta) deal(W_z_norm(md_hessian_analysis.evecs * (beta - beta_bar))^2 - covar_coeff * n, [], ...
                               2 * md_hessian_analysis.evecs' * z_prior_interface.Apply_W_z(md_hessian_analysis.evecs * (beta - beta_bar)), []);

        md_oed.Set_Covariance_Coefficient(covar_coeff);
        [beta_new, z_p] = md_oed.Generate_Seq_Optimal_Design_Con_v1(beta_0, alpha_d, betas, beta_bar, nonlcon);
        betas = [betas; beta_new];
        z_p = z_p(:, end);
        % disp(norm(z_p - z_bar) / norm(z_bar));
        % disp(beta_new - beta_bar);
        % disp(delta_beta);
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
    Jhat_DC_oed(p) = Jhat_hifi_fn(z_bar);
    fprintf('Objective of z_bar: \t%.3f\n', Jhat_DC_oed(p));
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

figure;
hold on;
plot(flip(z_lofi), "Color", [0.7 0.7 0.7], "DisplayName", "$\tilde{z}$", "LineWidth", 3);
plot(flip(z_hifi), "k-", "DisplayName", "$z_{\rm hifi}$", "LineWidth", 3);
plot(flip(z_bar), "b-", "DisplayName", "$\overline{z}$", "LineWidth", 3);
xlabel("$x$", "Interpreter", "latex");
ylabel("$c_0$", "Interpreter", "latex");
legend("location", "northeast", "Interpreter", "latex");

%%%%%%%%%%%%%%%%
% Evaluating Hi-Fi PDE at Hi-Fi Sol...

% Step 0:
% -------------
% Objective of z_lofi:    17.006
% Objective of z_hifi:    0.095
% Objective of z_proj:    0.134

% Step 1:
% -------------
% Objective of z_bar:     3.932
% Percent Improvement:    77.31%

% Step 2:
% -------------
% Objective of z_bar:     2.456
% Percent Improvement:    38.46%

% Step 3:
% -------------
% Objective of z_bar:     2.421
% Percent Improvement:    1.50%

% Step 4:
% -------------
% Objective of z_bar:     2.222
% Percent Improvement:    8.53%

% Step 5:
% -------------
% Objective of z_bar:     2.051
% Percent Improvement:    8.04%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
