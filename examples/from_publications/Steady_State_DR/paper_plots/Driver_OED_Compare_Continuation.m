% Get the OED Setup ready
addpath(genpath('..'));
addpath(genpath('../../../src'));
OED_Setup;

load("oed-results.mat");

% Perform Offline OED Computations - USES data_interface
md_oed = MD_OED(opt_prob_interface, data_interface, u_prior_interface, z_prior_interface, md_hessian_analysis);
md_oed.Offline_Computation();

%% Perform OED
N = 5;
Z = [];
D = [];
betas = [];
Jhat_noncnt = zeros(N, 1);
z_bars_noncnt = zeros(n, N);
beta_bars = zeros(num_evals, N);
beta_0 = randn(num_evals, 1);

for p = 1:N
    fprintf('\nStep %d:\n-------------\n', p);

    % Obtain Discrepancies
    Z = Z_oed(:, 1:p);
    D = D_oed(:, 1:p);
    data_interface.Set_Z_and_D(Z, D);

    % Perform Posterior Sampling (TODO: Reuse data)
    md_post_sampling = MD_Posterior_Sampling(data_interface, u_prior_interface, z_prior_interface);
    md_post_sampling.Compute_Posterior_Data(alpha_d, 1);

    % Obtain Optimal Solution Update via Continuation
    % num_continuation_steps = 1;
    md_cont_update = MD_Update(md_post_sampling, md_hessian_analysis);
    z_bar = md_cont_update.Posterior_Update_Mean();

    % Display Stats
    Jhat_noncnt(p) = opt_hifi.Jhat(z_bar);
    fprintf('Objective of z_bar: \t%.3f\n', Jhat_noncnt(p));

end

% Plot Objective Function over N
if true
    figure;
    hold on;
    xlim([0 N]);
    yline(Jhat_hifi, "k--", "DisplayName", "Hi-Fi", "LineWidth", 3, "Layer", "Bottom", "Alpha", 1);
    yline(Jhat_lofi, "r--", "DisplayName", "Lo-Fi", "LineWidth", 3, "Layer", "Bottom", "Alpha", 1);
    plot(0:N, [Jhat_lofi; Jhat_noncnt], ".-", "Color", "#EDB120", "DisplayName", "Post-optimality Linearization");
    plot(0:N, [Jhat_lofi; Jhat_DC_oed], ".-", "Color", "#77AC30", "DisplayName", "Continuation ($N_c = 3$)");
    xlabel("Evaluations ($N$)", "Interpreter", "latex");
    ylabel("High-fidelity objective", "Interpreter", "latex");
    legend("location", "east", "Interpreter", "latex");
    % title("Optimization Objective over Evals");
    saveas(gcf, 'ContinuationPlot', 'epsc');
end
