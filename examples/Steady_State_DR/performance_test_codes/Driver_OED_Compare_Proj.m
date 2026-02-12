% Get the OED Setup ready
addpath(genpath('..'));
addpath(genpath('../../../src'));
OED_Setup;

load("oed-results-con2-new.mat");
% pool = parpool("Threads", 8);

%% Perform Random Data Point Comparison
N = 5;
num_samples = 20;
Jhat_rand = zeros(N, num_samples);
rng(0);
Jhat_rand(1, :) = Jhat_DC_oed(1);

V = md_hessian_analysis.evecs;
% alpha_k_denom = trace(z_prior_interface.Apply_W_z_Inverse(z_prior_interface.M));
alpha_k_denom_proj = trace(V' * z_prior_interface.M * V) / W_z_norm(z_lofi)^2;

wb = waitbar(0, 'Starting Random');
for p = 2:N
    % fprintf('\nStep %d:\n-------------\n', p);
    % Generate random data points
    if p == 2
        alpha_k_num = M_z_norm(z_bars(:, p - 1) - z_lofi)^2;
    else
        alpha_k_num = M_z_norm(z_bars(:, p - 1) - z_bars(:, p - 2))^2;
    end

    alpha_k_proj = alpha_k_num / alpha_k_denom_proj;
    z_ps = z_bars(:, p - 1) + sqrt(alpha_k_proj) / W_z_norm(z_lofi) * V * randn(num_evals, num_samples);

    for i = 1:num_samples
        % disp(i)
        waitbar((num_samples * (p - 2) + i) / (num_samples * (N - 1)), wb, sprintf('Random Progress: %d %%', floor((num_samples * (p - 2) + i) / (num_samples * (N - 1)) * 100)));
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
    % waitbar((p-1) / (N-1), wb, sprintf('Random Progress: %d %%', floor((p-1) / (N-1) * 100)));
end
close(wb);

figure;
% yline(Jhat_hifi, "k--", "DisplayName", "Hi-Fi", "LineWidth", 3, "Layer", "Bottom", "Alpha", 1);
semilogy(1:N, Jhat_rand' - 3.086, ".", "Color", [1 1 1] - 0.5 * (1 - [0.7 0.8 0.9]), "MarkerSize", 25, "HandleVisibility", "off");
hold on;
semilogy(1:N, mean(Jhat_rand') - 3.086, ".--", "Color", [0.7 0.8 0.9], "MarkerSize", 25, "HandleVisibility", "off");
semilogy(0:N, [Jhat_lofi; Jhat_DC_oed] - 3.086, ".-", "Color", "#BAB86C", "DisplayName", "OED");
semilogy(1, Jhat_DC_oed(1) - 3.086, ".", "Color", [0.7 0.8 0.9], "MarkerSize", 25, "DisplayName", "Random");
yline(Jhat_lofi - 3.086, "r--", "DisplayName", "$\hat{J}(\tilde{z})$", "LineWidth", 3, "Layer", "Bottom", "Alpha", 1);
xlim([0 N]);
% plot(0:N, [Jhat_lofi; Jhat_DC_oed], ".-", "Color", "#BAB86C", "DisplayName", "DeltaCov OED");
xlabel("Evaluations ($N$)", "Interpreter", "latex");
ylabel("Objective $\hat{J}(\cdot)$", "Interpreter", "latex");
legend("location", "east", "Interpreter", "latex");
% title("Optimization Objective over Evals");
saveas(gcf, "contObjComparePlot.eps", "epsc");
