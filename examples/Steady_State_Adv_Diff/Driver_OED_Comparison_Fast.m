%% Set up
% clear;
% close all;
% clc;
addpath(genpath('../../src'));
load OED_Ensemble_Results.mat;
clear data_interface Z D;

p = 5;
num_post_samps = 1;
save_figures = false;

Jhat_lofi = opt_hifi.Jhat(z_lofi);
Jhat_hifi = opt_hifi.Jhat(z_hifi);

% Calculate Relative OED Error (Lambda Function for now)
oed_z_error_fn = @(z) sqrt((z - z_hifi)' * z_prior_interface.Apply_M_z(z - z_hifi)) / sqrt(z_hifi' * z_prior_interface.Apply_M_z(z_hifi));

z_error_0 = oed_z_error_fn(z_lofi);

Jhat_oed = zeros(p + 1, 1);
oed_z_error = zeros(p + 1, 1);

k = 1;
fprintf('\nStep %d:\n-------------\n', k);
Z = z_lofi;
D = Evaluate_Discrepancy(con_hifi, con_lofi, z_lofi);
data_interface = MD_Data_Interface_Diff(u_lofi, z_lofi);
data_interface.Set_Z_and_D(Z, D);
md_post_sampling = MD_Posterior_Sampling(data_interface, u_prior_interface, z_prior_interface);
md_post_sampling.Compute_Posterior_Data(alpha_d, num_post_samps);
md_update = MD_Update(md_post_sampling, md_hessian_analysis);
z_update_mean = md_update.Posterior_Update_Mean();
Jhat_oed(1) = opt_hifi.Jhat(z_update_mean);
oed_z_error(1) = oed_z_error_fn(z_update_mean);

tmp = Extract_mean_theta(z_prior_interface, md_post_sampling.post_data);
std_oed_mean_theta = zeros(length(tmp), p + 1);
std_oed_mean_theta(:, 1) = tmp;
std_oed_mean_z = zeros(length(z_update_mean), p + 1);
std_oed_mean_z(:, 1) = z_update_mean;
std_oed_Z = cell(p + 1, 1);
std_oed_Z{1} = z_lofi;

fprintf('Objective of z_bar: \t%.3f\n', Jhat_oed(k));
% fprintf('Rel. Err of z_bar: \t%.2f%%\n', 100 * oed_z_error(k));
% fprintf('Diff. w/ z_hifi obj.: \t%.2f%%\n', 100 * (Jhat_oed(k) - Jhat_hifi) / (Jhat_hifi));
fprintf('Percent Improvement: \t%.2f%%\n\n', 100 * (Jhat_lofi - Jhat_oed(k)) / (Jhat_lofi - Jhat_hifi));

for k = 2:p + 1
    fprintf('\nStep %d:\n-------------\n', k);

    data_interface = MD_Data_Interface_Diff(k - 1, 1, 'OED');
    data_interface.Load_Data();
    md_post_sampling = MD_Posterior_Sampling(data_interface, u_prior_interface, z_prior_interface);
    md_post_sampling.Compute_Posterior_Data(alpha_d, num_post_samps);
    md_update = MD_Update(md_post_sampling, md_hessian_analysis);
    z_update_mean = md_update.Posterior_Update_Mean();
    Jhat_oed(k) = opt_hifi.Jhat(z_update_mean);
    oed_z_error(k) = oed_z_error_fn(z_update_mean);

    std_oed_mean_theta(:, k) = Extract_mean_theta(z_prior_interface, md_post_sampling.post_data);
    std_oed_mean_z(:, k) = z_update_mean;
    std_oed_Z{k} = data_interface.Z;

    % Display Stats
    fprintf('Objective of z_bar: \t%.3f\n', Jhat_oed(k));
    % fprintf('Rel. Err of z_bar: \t%.2f%%\n', 100 * oed_z_error(k));
    % fprintf('Diff. w/ z_hifi obj.: \t%.2f%%\n', 100 * (Jhat_oed(k) - Jhat_hifi) / (Jhat_hifi));
    fprintf('Percent Improvement: \t%.2f%%\n\n', 100 * (Jhat_oed(k - 1) - Jhat_oed(k)) / (Jhat_oed(k - 1) - Jhat_hifi));

    if k == 20
        figure;
        hold on;
        plot(x, con_hifi.State_Solve(z_lofi), "r-", "DisplayName", "$S(\tilde{z})$");
        plot(x, con_hifi.State_Solve(z_update_mean), "g-", "DisplayName", "$S(\bar{z})$");
        plot(x, con_hifi.State_Solve(z_hifi), "k--", "DisplayName", "$S(z^*)$");
        % plot(x, obj.T, "k:", "DisplayName", "Target");
        % ylim(fixed_ylim);
        title("Optimization State (N = " + k + ")");
        legend("Location", "northwest", "interpreter", "latex");
        if save_figures
            saveas(gcf, "StdOED_N_" + k + ".svg");
        end
        figure;
        hold on;
        plot(x, z_lofi, "r-", "DisplayName", "Lo-Fi Sol.");
        plot(x, z_update_mean, "g-", "DisplayName", "Updated Sol.");
        plot(x, z_hifi, "k--", "DisplayName", "Hi-Fi Sol.");
        title("Updated Control (Iteration " + k + ")");
        legend("Location", "best");
        if save_figures
            saveas(gcf, "StdOED_Ctrl_N_" + k + ".svg");
        end
    end
end

old_oed = Jhat_oed;

figure;
hold on;
xlim([0 p]);
plot(xlim, 0 * xlim + Jhat_hifi, "k--", "DisplayName", "Hi-Fi");
plot(xlim, 0 * xlim + Jhat_lofi, "r--", "DisplayName", "Lo-Fi");
plot([0:p], [Jhat_lofi; old_oed(1:end - 1)], ".-", "Color", "#1F618D", "DisplayName", "Updated Sols");
xlabel("Evaluations ($N$)", "Interpreter", "latex");
ylabel("Objective $\hat{J}(\cdot)$", "Interpreter", "latex");
legend("location", "east", "Interpreter", "latex");
title("Optimization Objective over Evals");
if save_figures
    saveas(gcf, "Std_OED_Obj.svg");
end

save('Std_OED_Results.mat', 'std_oed_mean_theta', 'std_oed_mean_z', 'std_oed_Z');
