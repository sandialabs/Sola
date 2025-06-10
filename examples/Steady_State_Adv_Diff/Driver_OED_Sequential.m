% Clear Workspace and Add Interfaces to Path
% clear;
close all;
% clc;
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
obj = Adv_Diff_Objective(m, reg_coeff);
con_hifi = Adv_Diff_Constraint(m, diff_coeff, vel_coeff, robin_coeff);
con_lofi = Diff_Constraint(con_hifi);
opt_hifi = Reduced_Space_Optimization(obj, con_hifi);
opt_lofi = Reduced_Space_Optimization(obj, con_lofi);
x = con_lofi.x;

% Show initial objective
fprintf("\nStep 0:\n-------------\n");
Jhat_lofi = opt_hifi.Jhat(z_lofi);
Jhat_hifi = opt_hifi.Jhat(z_hifi);
fprintf('Objective of z_lofi: \t%.3f\n', Jhat_lofi);
fprintf('Objective of z_hifi: \t%.3f\n\n', Jhat_hifi);

% Set Data Interface (no data there yet, except for z_lofi/u_lofi)
data_interface = MD_Data_Interface_Diff(u_lofi, z_lofi);

% Generate Priors for u and z
alpha_u = (1 / 2)^2;
alpha_z = (1 / 100)^2;
alpha_d = 1.e-2;
u_prior_interface = MD_Elliptic_u_Prior_Interface_Diff(alpha_u, opt_lofi);
z_prior_interface = MD_Elliptic_z_Prior_Interface_Diff(alpha_z, opt_lofi);

% Error with z_hifi
oed_z_error_fn = @(z) sqrt((z - z_hifi)' * z_prior_interface.Apply_M_z(z - z_hifi)) / sqrt(z_hifi' * z_prior_interface.Apply_M_z(z_hifi));

% %% Hessian analysis
opt_prob_interface = MD_Opt_Prob_Interface_Sabl(opt_lofi, data_interface);
md_hessian_analysis = MD_Hessian_Analysis(opt_prob_interface, z_prior_interface);
num_evals = 4;
oversampling = 20;
md_hessian_analysis.Compute_Hessian_GEVP(data_interface.z_init, num_evals, oversampling);

% Perform Offline OED Computations - This is used to generate many random designs (to avoid OED in next steps)
Im = eye(m);
Mz = z_prior_interface.M;
B1 = @(x) opt_prob_interface.Apply_Solution_Operator_z_Jacobian_Transpose(opt_prob_interface.Apply_Misfit_Hessian([Im kron(Im, z_lofi' * Mz)] * x, u_lofi, z_lofi), z_lofi);
B2 = @(x) [zeros(m, m) kron(opt_prob_interface.Misfit_Gradient(u_lofi, z_lofi)', Mz)] * x;
B = @(x) B1(x) + B2(x);
PHinvB = @(x) md_hessian_analysis.Apply_Projected_RS_Hessian_Inverse(B(x));

reg_coeff = 1.e-6;
beta_0 = randn(num_evals, 1);
oed_interface = MD_OED_Interface_Diff(data_interface, con_lofi);

% Plot low-fidelity and high-fidelity states
if show_figures
    figure;
    hold on;
    plot(x, con_hifi.State_Solve(z_lofi), "r-", "DisplayName", "$S(\tilde{z})$");
    plot(x, con_hifi.State_Solve(z_hifi), "k--", "DisplayName", "$S(z^*)$");
    % plot(x, obj.T, "k:", "DisplayName", "Target");
    title("Seq-OED State (Iteration 0)");
    fixed_ylim = ylim;
    legend("Location", "northwest", "interpreter", "latex");
    if save_figures
        saveas(gcf, "SeqOED_N_0.png");
    end
end

%% Iterate for each data point
N = 6;
Jhat_oed = zeros(N, 1);
oed_z_error = zeros(N, 1);
Z = [];
D = [];
betas = [];
z_bar = z_lofi;

m = length(z_lofi);
seq_oed_mean_theta = zeros(m * (m + 1), N);
seq_oed_mean_z = zeros(m, N);
seq_oed_Z = cell(N, 1);

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
        [beta_new, z_p] = md_oed.Generate_Seq_Optimal_Design(beta_0, alpha_d, reg_coeff, betas);
        betas = [betas; beta_new];
        z_p = z_p(:, end);
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

    % Obtain Optimal Solution Update
    md_update = MD_Update(md_post_sampling, md_hessian_analysis);
    z_bar = md_update.Posterior_Update_Mean();
    theta_post = Extract_mean_theta(md_post_sampling.post_data);
    % z_bar = z_lofi - PHinvB(theta_post);

    % Display Stats
    Jhat_oed(p) = opt_hifi.Jhat(z_bar);
    oed_z_error(p) = oed_z_error_fn(z_bar);
    fprintf('Objective of z_bar: \t%.3f\n', Jhat_oed(p));
    % fprintf('Rel. Err of z_bar: \t%.2f%%\n', 100 * oed_z_error(p));
    % fprintf('Diff. w/ z_hifi obj.: \t%.2f%%\n', 100 * (Jhat_oed(p) - Jhat_hifi) / (Jhat_hifi));
    if p == 1
        fprintf('Percent Improvement: \t%.2f%%\n\n', 100 * (Jhat_lofi - Jhat_oed(p)) / (Jhat_lofi - Jhat_hifi));
    else
        fprintf('Percent Improvement: \t%.2f%%\n\n', 100 * (Jhat_oed(p - 1) - Jhat_oed(p)) / (Jhat_oed(p - 1) - Jhat_hifi));
    end

    % Plot Low-Fidelity, High-fidelity and Updated States
    if show_figures
        figure;
        hold on;
        plot(x, con_hifi.State_Solve(z_lofi), "r-", "DisplayName", "$S(\tilde{z})$");
        plot(x, con_hifi.State_Solve(z_bar), "b-", "DisplayName", "$S(\bar{z})$");
        plot(x, con_hifi.State_Solve(z_hifi), "k--", "DisplayName", "$S(z^*)$");
        % plot(x, obj.T, "k:", "DisplayName", "Target");
        ylim(fixed_ylim);
        title("Optimization State (N = " + p + ")");
        legend("Location", "northwest", "interpreter", "latex");
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

    seq_oed_mean_theta(:, p) = theta_post;
    seq_oed_mean_z(:, p) = z_bar;

end

% Plot Objective Function over N
show_figures = true;
if show_figures
    figure;
    hold on;
    xlim([0 N]);
    yline(Jhat_hifi, "k--", "DisplayName", "Hi-Fi", "LineWidth", 3, "Layer", "Bottom", "Alpha", 1);
    yline(Jhat_lofi, "r--", "DisplayName", "Lo-Fi", "LineWidth", 3, "Layer", "Bottom", "Alpha", 1);
    try
        yline(Jhat_HDSA, "b--", "DisplayName", "Best-HDSA", "LineWidth", 3, "Layer", "Bottom", "Alpha", 1);
        plot(0:N, [Jhat_lofi; old_oed(1:N)], ".-", "Color", "#1F618D", "DisplayName", "Standard OED");
    catch ME
        if ~strcmp(ME.identifier, 'MATLAB:UndefinedFunction')
            rethrow(ME);
        end
    end
    plot(0:N, [Jhat_lofi; Jhat_oed], ".-", "Color", "#00C83A", "DisplayName", "Sequential OED");
    xlabel("Evaluations ($N$)", "Interpreter", "latex");
    ylabel("Objective $\hat{J}(\cdot)$", "Interpreter", "latex");
    legend("location", "east", "Interpreter", "latex");
    title("Optimization Objective over Evals");
    if save_figures
        saveas(gcf, "SeqOED_Objs.png");
    end
end

save('Seq_OED_Results.mat', 'seq_oed_mean_theta', 'seq_oed_mean_z', 'seq_oed_Z');
