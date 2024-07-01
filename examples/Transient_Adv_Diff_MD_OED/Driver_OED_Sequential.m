% Clear Workspace and Add Interfaces to Path
addpath(genpath('../../src'));
rng(0);

% Set Default Font Axes and Line Width
set(0, "DefaultAxesFontSize", 20);
set(0, "DefaultLineLineWidth", 3);
set(0, "DefaultLineMarkerSize", 20);

show_figures = false;
save_figures = false;

% Retrieve Model Parameters (D, Z, diff/reg/react_coeff, m, u_lofi, z_hifi/lofi; remove Z and D though)
load Optimization_Results.mat;
clear Z D;

% Set Hi-Fi and Lo-Fi Objectives and Constraints
n_z = num_space_control_nodes * (n_t - 1);
obj = Adv_Diff_Gaussian_Source_Objective(n_y, n_z, T, n_t, num_space_control_nodes, reg_coeff);
con_lofi = Adv_Diff_Gaussian_Source_Constraint(n_y, n_z, T, n_t, num_space_control_nodes, diff_coeff, vel_coeff_lofi);
opt_lofi = Reduced_Space_Optimization(obj, con_lofi);
con_hifi = Adv_Diff_Gaussian_Source_Constraint(n_y, n_z, T, n_t, num_space_control_nodes, diff_coeff, vel_coeff_hifi);
opt_hifi = Reduced_Space_Optimization(obj, con_hifi);
x = con_lofi.x;
t = con_lofi.t_mesh;
terminal_target = obj.Evaluate_Target(T, x);
index_terminal = @(full_u) full_u((n_t - 1) * n_y + 1:end, :);

% Show initial objective
fprintf("\nStep 0:\n-------------");
Jhat_lofi = opt_hifi.Jhat(z_lofi);
Jhat_hifi = opt_hifi.Jhat(z_hifi);
fprintf('Objective of z_lofi: \t%.6f\n', Jhat_lofi);
fprintf('Objective of z_hifi: \t%.6f\n\n', Jhat_hifi);

% Set Data Interface (no data there yet, except for z_lofi/u_lofi)
data_interface = MD_Data_Interface_Adv_Diff(u_lofi, z_lofi);

% Set Transient Prior
beta_t = 50;
beta_i = 1.e5;
transient_prior_cov = MD_Transient_Prior_Covariance_Sabl(beta_t, beta_i, T, n_t, n_y);

% Generate Priors for u and z
alpha_u = (1 / 1)^2;
alpha_d = 1.e-6;
u_prior_interface = MD_Transient_Elliptic_u_Prior_Interface_Adv_Diff(alpha_u, transient_prior_cov, opt_lofi);
z_prior_interface = MD_z_Prior_Interface_Adv_Diff(obj);

% Error with z_hifi (normalized w/ z_lofi)
normalization = Control_Norm(z_hifi, con_lofi);
oed_z_error_fn = @(z) Control_Norm(z - z_hifi, con_lofi) / normalization;

% %% Hessian analysis
opt_prob_interface = MD_Opt_Prob_Interface_Sabl(opt_lofi, data_interface);
md_hessian_analysis = MD_Hessian_Analysis(opt_prob_interface, z_prior_interface);
num_evals = 4; % CHANGED FROM 26
oversampling = 20;
md_hessian_analysis.Compute_Hessian_GEVP(data_interface.z_init, num_evals, oversampling);

% Perform Offline OED Computations - This is used to generate many random designs (to avoid OED in next steps)
reg_coeff = 1.e-5;
beta_0 = randn(num_evals, 1);
oed_interface = MD_OED_Interface_Adv_Diff(data_interface, obj);

% Matrix Projection operator for sequential OED
project_to_Z = @(Z, z_oed) Z * (Z \ z_oed);

% Plot low-fidelity and high-fidelity states

if show_figures
    figure;
    hold on;
    plot(x, terminal_target, "c-", "DisplayName", "Target");
    plot(x, index_terminal(con_hifi.State_Solve(z_lofi)), "r-", "DisplayName", "$S(\tilde{z})$");
    plot(x, index_terminal(con_hifi.State_Solve(z_hifi)), "k--", "DisplayName", "$S(z^*)$");
    title("Seq-OED Terminal State (Iteration 0)");
    fixed_ylim = 2 * ylim - mean(ylim);
    ylim(fixed_ylim);
    legend("Location", "southeast", "interpreter", "latex");
    if save_figures
        saveas(gcf, "SeqOED_N_0.png");
    end
end

%% Iterate for each data point
N = 10;
Jhat_oed = zeros(N, 1);
oed_z_error = zeros(N, 1);
Z = [];
D = [];
betas = [];
z_bar = z_lofi;

for p = 1:N
    % Update Data Interface (with prior center)
    fprintf('\nStep %d:\n-------------\n', p);
    data_interface.Update_z_opt(z_bar);

    % Sequential OED
    md_oed = MD_OED_Seq(opt_prob_interface, data_interface, u_prior_interface, z_prior_interface, md_hessian_analysis, oed_interface);
    md_oed.Offline_Computation();

    % Set Parameters for OED
    if p == 1
        z_p = z_lofi;
        % betas = [betas; 0*beta_0]; % This is for the updated sequential OED
    else
        [beta_new, z_p] = md_oed.Generate_Seq_Optimal_Design(beta_0, alpha_d, reg_coeff, betas);
        betas = [betas; beta_new];
        z_p = z_p(:, end); % Redundancy for standard OED.
    end

    % Obtain Discrepancies
    Z = [Z z_p];
    D_p = Evaluate_Discrepancy(con_hifi, con_lofi, z_p);
    D = [D D_p];
    data_interface.Set_Z_and_D(Z, D);

    % New Idea: revert z_opt (prior term) to span of Z
    % z_h = project_to_Z(Z, z_bar);
    % disp(norm(z_bar - z_h))
    % data_interface.Update_z_opt(z_h);

    % Perform Posterior Sampling
    md_post_sampling = MD_Posterior_Sampling(data_interface, u_prior_interface, z_prior_interface);
    md_post_sampling.Compute_Posterior_Data(alpha_d, 1);

    % Obtain Optimal Solution Update
    md_update = MD_Update(md_post_sampling, md_hessian_analysis);
    z_bar = md_update.Posterior_Update_Mean();

    % Display Stats
    Jhat_oed(p) = opt_hifi.Jhat(z_bar);
    oed_z_error(p) = oed_z_error_fn(z_bar);
    fprintf('Objective of z_bar: \t%.5f\n', Jhat_oed(p));
    % fprintf('Rel. Err w/ z_hifi: \t%.3f\n', oed_z_error(p));
    fprintf('Diff. w/ z_hifi obj.: \t%.2f%%\n', 100 * (Jhat_oed(p) - Jhat_hifi) / (Jhat_lofi - Jhat_hifi));
    if p == 1
        fprintf('Percent Improvement: \t%.2f%%\n\n', 100 * (Jhat_lofi - Jhat_oed(p)) / (Jhat_lofi - Jhat_hifi));
    else
        fprintf('Percent Improvement: \t%.2f%%\n\n', 100 * (Jhat_oed(p - 1) - Jhat_oed(p)) / (Jhat_oed(p - 1) - Jhat_hifi));
    end

    % Plot Low-Fidelity, High-fidelity and Updated States
    if show_figures
        figure;
        hold on;
        plot(x, index_terminal(con_hifi.State_Solve(z_lofi)), "r-", "DisplayName", "$S(\tilde{z})$");
        plot(x, index_terminal(con_hifi.State_Solve(z_hifi)), "k--", "DisplayName", "$S(z^*)$");
        plot(x, index_terminal(con_hifi.State_Solve(z_bar)), "b-", "DisplayName", "$S(\bar{z})$");
        ylim(fixed_ylim);
        title("Seq-OED Terminal State (Iteration " + p + ")");
        legend("Location", "southeast", "interpreter", "latex");
        if save_figures
            saveas(gcf, "SeqOED_N_" + p + ".png");
        end
        % figure;
        % hold on;
        % plot(x, z_lofi, "r-", "DisplayName", "Lo-Fi Sol.");
        % plot(x, z_hifi, "k--", "DisplayName", "Hi-Fi Sol.");
        % plot(x, z_bar, "b-", "DisplayName", "Updated Sol.");
        % title("Lo-Fi & Hi-Fi Controls (Iteration " + p + ")");
        % legend("Location", "best");
    end

end

% Plot Objective Function over N
show_figures = true;
if show_figures
    figure;
    hold on;
    xlim([0 N]);
    yline(Jhat_hifi, "k--", "DisplayName", "Hi-Fi", "LineWidth", 3, "Layer", "Bottom", "Alpha", 1);
    yline(Jhat_lofi, "r--", "DisplayName", "Lo-Fi", "LineWidth", 3, "Layer", "Bottom", "Alpha", 1);
    % plot(0:N, [Jhat_lofi; old_oed(1:N)], ".-", "Color", "#1F618D", "DisplayName", "Standard OED")
    plot(0:N, [Jhat_lofi; Jhat_oed], ".-", "Color", "#00C83A", "DisplayName", "Sequential OED");
    xlabel("Evaluations ($N$)", "Interpreter", "latex");
    ylabel("Objective $\hat{J}(\cdot)$", "Interpreter", "latex");
    legend("location", "east", "Interpreter", "latex");
    title("Optimization Objective over Evals");

    figure;
    hold on;
    xlim([0 N]);
    yline(oed_z_error_fn(z_lofi), "r--", "DisplayName", "Lo-Fi", "LineWidth", 3, "Layer", "Bottom", "Alpha", 1);
    plot(0:N, [oed_z_error_fn(z_lofi); oed_z_error], ".-", "Color", "#00C83A", "DisplayName", "Sequential OED");
    xlabel("Evaluations ($N$)", "Interpreter", "latex");
    ylabel("Rel. Err w/ $\bar{z}$", "Interpreter", "latex");
    legend("location", "northwest", "Interpreter", "latex");
    title("Relative Error of Sol. over Evals");

    figure;
    u_full_lofi = reshape(con_hifi.State_Solve(z_lofi), n_y, n_t);
    u_full_hifi = reshape(con_hifi.State_Solve(z_hifi), n_y, n_t);
    u_full_bar = reshape(con_hifi.State_Solve(z_bar), n_y, n_t);
    fixed_ylim_anim = [min([u_full_lofi u_full_hifi u_full_bar], [], "all") max([u_full_lofi u_full_hifi u_full_bar], [], "all")];
    for k = 1:n_t
        plot(x, u_full_lofi(:, k), "r-", "DisplayName", "$S(\tilde{z})$");
        hold on;
        plot(x, u_full_hifi(:, k), "k--", "DisplayName", "$S(z^*)$");
        plot(x, u_full_bar(:, k), "b-", "DisplayName", "$S(\bar{z})$");
        ylim(fixed_ylim_anim);
        legend("Location", "east", "interpreter", "latex");
        hold off;
        pause(.05);
    end
end
