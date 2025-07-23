% Get the OED Setup ready
addpath(genpath('..'));
addpath(genpath('../../../src'));
OED_Setup;

% Perform Offline OED Computations - USES data_interface
md_oed = MD_OED_DeltaCov(opt_prob_interface, data_interface, u_prior_interface, z_prior_interface, md_hessian_analysis, oed_interface);
md_oed.Offline_Computation();

%% Perform OED
N = 5;
Z_oed = [];
D_oed = [];
betas_oed = [];
Jhat_DC_oed = zeros(N, 1);
oed_reg_coeff = 1.e-2;
z_bars = zeros(n, N);
beta_0 = randn(num_evals, 1);

wb = waitbar(0, 'Starting OED');
for p = 1:N
    waitbar(p / N, wb, sprintf('OED Progress: %d %%', floor(p / N * 100)));
    % fprintf('\nStep %d:\n-------------\n', p);

    % Set Parameters for OED
    if p == 1
        z_p = z_lofi;
    else
        if p == 2
            covar_coeff = W_z_norm(z_bar - z_lofi)^2 / n;
        else
            covar_coeff = W_z_norm(z_bar - z_bars(:, p - 2))^2 / n;
        end
        md_oed.Set_Covariance_Coefficient(covar_coeff);
        beta_0 = betas_cont(:, end);
        [beta_new, z_p] = md_oed.Generate_Seq_Optimal_Design(beta_0, alpha_d, oed_reg_coeff, betas_oed, betas_cont(:, end));
        betas_oed = [betas_oed; beta_new];
        z_p = z_p(:, end);
        % disp(norm(z_p - z_bar) / norm(z_bar));
    end

    % Obtain Discrepancies
    Z_oed = [Z_oed z_p];
    D_p = Evaluate_Discrepancy(con_hifi, con_lofi, z_p);
    D_oed = [D_oed D_p];
    data_interface.Set_Z_and_D(Z_oed, D_oed);

    % Perform Posterior Sampling (TODO: Reuse data)
    md_post_sampling = MD_Posterior_Sampling(data_interface, u_prior_interface, z_prior_interface);
    md_post_sampling.Compute_Posterior_Data(alpha_d, 1);

    % Obtain Optimal Solution Update via Continuation
    num_continuation_steps = 3;
    md_cont_update = MD_Continuation_Update(md_post_sampling, md_hessian_analysis, num_continuation_steps);
    [u_cont, z_cont, betas_cont] = md_cont_update.Posterior_Update_Mean_PC_beta();
    z_bar = z_cont(:, end);
    z_bars(:, p) = z_bar;
    Jhat_DC_oed(p) = opt_hifi.Jhat(z_bar);
end
close(wb);

%% Perform Random Data Point Comparison
N = 5;
num_samples = 20;
Jhat_rand = zeros(N, num_samples);
rng(0);
Jhat_rand(1, :) = Jhat_DC_oed(1);

wb = waitbar(0, 'Starting Random');
for p = 2:N
    % fprintf('\nStep %d:\n-------------\n', p);
    % Generate random data points
    if p == 2
        covar_coeff = W_z_norm(z_bars(:, p - 1) - z_lofi)^2 / n;
    else
        covar_coeff = W_z_norm(z_bars(:, p - 1) - z_bars(:, p - 2))^2 / n;
    end
    z_ps = z_bar + sqrt(covar_coeff) * z_prior_interface.Sample_with_Covariance_W_z_Inverse(num_samples);

    for i = 1:num_samples
        % disp(i)
        waitbar(i / num_samples * (p - 1) / (N - 1), wb, sprintf('Random Progress: %d %%', floor(i / num_samples * (p - 1) / (N - 1) * 100)));
        % Obtain random data point
        z_p = z_ps(:, i);
        D_p = Evaluate_Discrepancy(con_hifi, con_lofi, z_p);

        % Obtain Discrepancies
        Z_rand = [Z_oed(:, 1:p - 1) z_p];
        D_rand = [D_oed(:, 1:p - 1) D_p];
        data_interface.Set_Z_and_D(Z_rand, D_rand);

        % Perform Posterior Sampling (TODO: Reuse data)
        md_post_sampling = MD_Posterior_Sampling(data_interface, u_prior_interface, z_prior_interface);
        md_post_sampling.Compute_Posterior_Data(alpha_d, 1);

        % Obtain Optimal Solution Update via Continuation
        num_continuation_steps = 3;
        md_cont_update = MD_Continuation_Update(md_post_sampling, md_hessian_analysis, num_continuation_steps);
        [u_cont, z_cont, betas_cont] = md_cont_update.Posterior_Update_Mean_PC_beta();
        Jhat_rand(p, i) = opt_hifi.Jhat(z_cont(:, end));
    end
end
close(wb);

figure;
% yline(Jhat_hifi, "k--", "DisplayName", "Hi-Fi", "LineWidth", 3, "Layer", "Bottom", "Alpha", 1);
semilogy(1:N, Jhat_rand' - 3.086, ".", "Color", [0.7 0.8 0.9], "MarkerSize", 25, "HandleVisibility", "off");
hold on;
semilogy(0:N, [Jhat_lofi; Jhat_DC_oed] - 3.086, ".-", "Color", "#BAB86C", "DisplayName", "OED");
semilogy(1, Jhat_DC_oed(1) - 3.086, ".", "Color", [0.7 0.8 0.9], "MarkerSize", 25, "DisplayName", "Random");
yline(Jhat_lofi - 3.086, "r--", "DisplayName", "$\hat{J}(\tilde{z})$", "LineWidth", 3, "Layer", "Bottom", "Alpha", 1);
xlim([0 N]);
% plot(0:N, [Jhat_lofi; Jhat_DC_oed], ".-", "Color", "#BAB86C", "DisplayName", "DeltaCov OED");
xlabel("Evaluations ($N$)", "Interpreter", "latex");
ylabel("Objective $\hat{J}(\cdot)$", "Interpreter", "latex");
legend("location", "east", "Interpreter", "latex");
title("Optimization Objective over Evals");
saveas(gcf, "oed-rand-2.eps", "epsc");
