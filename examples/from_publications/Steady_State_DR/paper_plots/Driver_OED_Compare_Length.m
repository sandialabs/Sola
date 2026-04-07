% Get the OED Setup ready
addpath(genpath('..'));
addpath(genpath('../../../src'));
OED_Setup;

% Perform Offline OED Computations - USES data_interface
md_oed = MD_OED(opt_prob_interface, data_interface, u_prior_interface, z_prior_interface, md_hessian_analysis);
md_oed.Offline_Computation();
alpha_k_denom = trace(z_prior_interface.Apply_W_z_Inverse(z_prior_interface.M));

rng(0);

%% Perform Random Data Point Comparison (p = 1)
N = 6;
p_seq = 1;
Z = [];
D = [];
betas = [];
Jhat_DC_oed_1 = zeros(N, 1);
z_bars = zeros(n, N);
beta_bars = zeros(num_evals, N);
beta_0 = randn(num_evals, p_seq);

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
    end

    % Obtain Discrepancies
    Z = [Z z_p];
    D_p = Evaluate_Discrepancy(con_hifi, con_lofi, z_p);
    D = [D D_p];
    data_interface.Set_Z_and_D(Z, D);

    % Perform Posterior Sampling
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
    Jhat_DC_oed_1(p) = opt_hifi.Jhat(z_bar);
    fprintf('Objective of z_bar: \t%.4f\n', Jhat_DC_oed_1(p));

end

%% Perform Random Data Point Comparison (p = 2)
N = 6;
p_seq = 2;
Z = [Z(:, 1)];
D = [D(:, 1)];
betas = [];
Jhat_DC_oed_2 = zeros(N, 1);
z_bars = zeros(n, N);
beta_bars = zeros(num_evals, N);
beta_0 = randn(num_evals * p_seq, 1);

for p = 1:(N / p_seq)
    fprintf('\nStep %d:\n-------------\n', 2 + (p - 1) * (p_seq));
    if p == 1
        alpha_k = M_z_norm(z_lofi)^2 / alpha_k_denom;
        md_oed.Set_Covariance_Coefficient(alpha_k);
        [beta_new, z_p] = md_oed.Generate_Seq_Optimal_Design(beta_0(num_evals + 1:end), alpha_d, betas, zeros(num_evals, 1), alpha_k_num);
        betas = [betas; beta_new];
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
    end

    % Obtain Discrepancies
    Z = [Z z_p];
    D_p = Evaluate_Discrepancy(con_hifi, con_lofi, z_p);
    D = [D D_p];
    data_interface.Set_Z_and_D(Z, D);

    % Perform Posterior Sampling
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
    Jhat_DC_oed_2(p) = opt_hifi.Jhat(z_bar);
    fprintf('Objective of z_bar: \t%.4f\n', Jhat_DC_oed_2(p));

end

%% Perform Random Data Point Comparison (p = 3)
N = 6;
p_seq = 3;
Z = [Z(:, 1)];
D = [D(:, 1)];
betas = [];
Jhat_DC_oed_3 = zeros(N, 1);
z_bars = zeros(n, N);
beta_bars = zeros(num_evals, N);
beta_0 = randn(num_evals * p_seq, 1);

for p = 1:(N / p_seq)
    fprintf('\nStep %d:\n-------------\n', 2 + (p - 1) * (p_seq));
    if p == 1
        alpha_k = M_z_norm(z_lofi)^2 / alpha_k_denom;
        md_oed.Set_Covariance_Coefficient(alpha_k);
        [beta_new, z_p] = md_oed.Generate_Seq_Optimal_Design(beta_0(num_evals + 1:end), alpha_d, betas, zeros(num_evals, 1), alpha_k_num);
        betas = [betas; beta_new];
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
    end

    % Obtain Discrepancies
    Z = [Z z_p];
    D_p = Evaluate_Discrepancy(con_hifi, con_lofi, z_p);
    D = [D D_p];
    data_interface.Set_Z_and_D(Z, D);

    % Perform Posterior Sampling
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
    Jhat_DC_oed_3(p) = opt_hifi.Jhat(z_bar);
    fprintf('Objective of z_bar: \t%.4f\n', Jhat_DC_oed_2(p));

end

%% Perform Random Data Point Comparison (p = 6)
N = 6;
p_seq = 6;
Z = [Z(:, 1)];
D = [D(:, 1)];
betas = [];
Jhat_DC_oed_6 = zeros(N, 1);
z_bars = zeros(n, N);
beta_bars = zeros(num_evals, N);
beta_0 = randn(num_evals * p_seq, 1);

for p = 1:(N / p_seq)
    fprintf('\nStep %d:\n-------------\n', 2 + (p - 1) * (p_seq));
    if p == 1
        alpha_k = M_z_norm(z_lofi)^2 / alpha_k_denom;
        md_oed.Set_Covariance_Coefficient(alpha_k);
        [beta_new, z_p] = md_oed.Generate_Seq_Optimal_Design(beta_0(num_evals + 1:end), alpha_d, betas, zeros(num_evals, 1), alpha_k_num);
        betas = [betas; beta_new];
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
    end

    % Obtain Discrepancies
    Z = [Z z_p];
    D_p = Evaluate_Discrepancy(con_hifi, con_lofi, z_p);
    D = [D D_p];
    data_interface.Set_Z_and_D(Z, D);

    % Perform Posterior Sampling
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
    Jhat_DC_oed_6(p) = opt_hifi.Jhat(z_bar);
    fprintf('Objective of z_bar: \t%.4f\n', Jhat_DC_oed_2(p));

end

figure;
semilogy(0:N, [Jhat_lofi; Jhat_DC_oed_1], ".-", "Color", "#BAB86C", "DisplayName", "$p=1$");
hold on;
plot([0 2:2:N], [Jhat_lofi; nonzeros(Jhat_DC_oed_2)], "r.-", "DisplayName", "$p=2$");
plot([0 3:3:N], [Jhat_lofi; nonzeros(Jhat_DC_oed_3)], "g.-", "DisplayName", "$p=3$");
plot([0 N], [Jhat_lofi; nonzeros(Jhat_DC_oed_6)], "k.-", "DisplayName", "$p=6$");
yline(Jhat_lofi, "r--", "DisplayName", "$\hat{J}(\tilde{z})$", "LineWidth", 3, "Layer", "Bottom", "Alpha", 1);
xlim([0 N]);
ylim([3.0855 inf]);
xlabel("Evaluations ($N$)", "Interpreter", "latex");
ylabel("High-fidelity objective", "Interpreter", "latex");
legend("location", "best", "Interpreter", "latex");
