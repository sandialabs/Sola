%% Set up
% clear;
% close all;
% clc;
addpath(genpath('../../src'));
load OED_Ensemble_Results.mat;
clear data_interface Z D;

% Set Default Font Axes and Line Width
set(0, "DefaultAxesFontSize", 20);
set(0, "DefaultLineLineWidth", 3);
set(0, "DefaultLineMarkerSize", 20);

N = 5;
num_post_samps = 1;
save_figures = false;

% Print Diagnostics
Jhat_lofi = opt_hifi.Jhat(z_lofi);
Jhat_hifi = opt_hifi.Jhat(z_hifi);
fprintf('\nStep %d:\n-------------\n', 0);
fprintf("\nObjective at z_lofi: \t" + Jhat_lofi);
fprintf("\nObjective at z_hifi: \t" + Jhat_hifi);
fprintf("\n\n");

% Calculate Relative OED Error (Lambda Function for now)
oed_z_error_fn = @(z) sqrt((z - z_hifi)' * z_prior_interface.Apply_M_z(z - z_hifi)) / sqrt(z_hifi' * z_prior_interface.Apply_M_z(z_hifi));
z_error_0 = oed_z_error_fn(z_lofi);
Jhat_std_oed = zeros(N, 1);
oed_z_error = zeros(N, 1);

% Saving to File
std_oed_mean_theta = cell(N, 1);
std_oed_mean_z = cell(N, 1);
std_oed_Z = cell(N, 1);

oedMethod = "live-stdOED"; % options: load-stdOED, live-stdOED, seqOED

% Some prep for other OED methods
seq_oed_Z = load("Seq_OED_Results.mat", "seq_oed_Z").seq_oed_Z;
Z = [];
D = [];
betas = [];
data_interface = MD_Data_Interface_Diff(u_lofi, z_lofi);
beta_0 = randn(num_evals, 1);

for k = 1:N
    % Run step k
    fprintf('\nStep %d:\n-------------\n', k);

    switch oedMethod
        case "load-stdOED"
            data_interface = MD_Data_Interface_Diff(k, 1, 'OED');
            data_interface.Load_Data();
        case "seqOED"
            data_interface.Set_Z_and_D(seq_oed_Z{k}, Evaluate_Discrepancy(con_hifi, con_lofi, seq_oed_Z{k}));
        case "live-stdOED"
            % Set Parameters for OED
            if k == 1
                z_p = z_lofi;
            else
                [beta_new, z_p] = md_oed.Generate_Seq_Optimal_Design(beta_0, alpha_d, oed_reg_coeff, betas);
                betas = [betas; beta_new];
                z_p = z_p(:, end);
            end

            % Obtain Discrepancies
            Z = [Z z_p];
            D_p = Evaluate_Discrepancy(con_hifi, con_lofi, z_p);
            D = [D D_p];
            data_interface.Set_Z_and_D(Z, D);
    end

    md_post_sampling = MD_Posterior_Sampling(data_interface, u_prior_interface, z_prior_interface);
    md_post_sampling.Compute_Posterior_Data(alpha_d, num_post_samps);
    % md_update = MD_Update(md_post_sampling, md_hessian_analysis);
    % z_update_mean = md_update.Posterior_Update_Mean();
    num_continuation_steps = 1;
    md_cont_update = MD_Continuation_Update(md_post_sampling, md_hessian_analysis, num_continuation_steps);
    [~, z_cont, ~] = md_cont_update.Posterior_Update_Mean_PC_beta();
    z_update_mean = z_cont(:, end);
    Jhat_std_oed(k) = opt_hifi.Jhat(z_update_mean);
    oed_z_error(k) = oed_z_error_fn(z_update_mean);

    % Saving to File
    std_oed_mean_theta{k} = Extract_mean_theta(md_post_sampling.post_data);
    std_oed_mean_z{k} = z_update_mean;
    std_oed_Z{k} = data_interface.Z;

    % Display Stats
    fprintf('Objective of z_bar: \t%.3f\n', Jhat_std_oed(k));
    if k > 1
        fprintf('Percent Improvement: \t%.2f%%\n\n', 100 * (Jhat_std_oed(k - 1) - Jhat_std_oed(k)) / (Jhat_std_oed(k - 1) - Jhat_hifi));
    end
end

figure;
hold on;
xlim([0 N]);
plot(xlim, 0 * xlim + Jhat_hifi, "k--", "DisplayName", "Hi-Fi");
plot(xlim, 0 * xlim + Jhat_lofi, "r--", "DisplayName", "Lo-Fi");
plot([0:N], [Jhat_lofi; Jhat_std_oed(1:end)], ".-", "Color", "#1F618D", "DisplayName", "Updated Sols");
xlabel("Evaluations ($N$)", "Interpreter", "latex");
ylabel("Objective $\hat{J}(\cdot)$", "Interpreter", "latex");
legend("location", "east", "Interpreter", "latex");
title("Optimization Objective over Evals");
if save_figures
    saveas(gcf, "Std_OED_Obj.svg");
end

save('Std_OED_Results.mat', 'Jhat_std_oed', 'std_oed_mean_theta', 'std_oed_mean_z', 'std_oed_Z');
