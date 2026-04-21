%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
close all;
clc;

load HDSA_Results.mat;

figure;
hold on;
plot(t, u_lofi_t(:, 1), 'color', 'black', 'LineWidth', 3);
plot(t, u_hifi_t(:, 1), 'color', 'cyan', 'LineWidth', 3);
plot(t, u_update_t(:, 1), '--', 'color', 'red', 'LineWidth', 3);
legend({'Low-fidelity controlled state', 'High-fidelity controlled state' 'Update controlled state'}, 'Location', 'northwest');
set(gca, 'fontsize', 24);
exportgraphics(gcf, 'States.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none');

figure;
hold on;
plot(t(2:end), z_tilde, 'color', 'black', 'LineWidth', 3);
plot(t(2:end), z_star, 'color', 'cyan', 'LineWidth', 3);
plot(t(2:end), z_update_mean, '--', 'color', 'red', 'LineWidth', 3);
for k = 1:num_post_samples
    plot(t(2:end), z_update_samples(:, k), 'color', [.9, .9, .9], 'LineWidth', 3);
end
plot(t(2:end), z_tilde, 'color', 'black', 'LineWidth', 3);
plot(t(2:end), z_star, 'color', 'cyan', 'LineWidth', 3);
plot(t(2:end), z_update_mean, '--', 'color', 'red', 'LineWidth', 3);
legend({'Low-fidelity control', 'High-fidelity control', 'Update'}, 'Location', 'northwest');
set(gca, 'fontsize', 24);
exportgraphics(gcf, 'Controls.pdf', 'ContentType', 'vector', 'BackgroundColor', 'none');
