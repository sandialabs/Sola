clear;
close all;
clc;

addpath(genpath('../../src'));
load('HDSA_Results.mat');

working_path = pwd;
write_path = '~/Desktop/junk'; % '~/Desktop/Model_Discrepancy_Sampling/figures/diff_react';

figure;
hold on;
plot(x, z_lofi, 'LineWidth', 3, 'color', 'black');
xlabel('$x$', 'Interpreter', 'latex');
ylabel('Source');
legend('$\tilde{z}$', 'Interpreter', 'latex', 'FontSize', 18, 'Location', 'northwest');
set(gca, 'fontsize', 18);
cd(write_path);
saveas(gca, 'lofi_source', 'epsc');
cd(working_path);

figure;
hold on;
plot(x, obj.T, '--', 'LineWidth', 3, 'color', 'magenta');
plot(x, con_lofi.State_Solve(z_lofi), 'LineWidth', 3, 'color', 'green');
plot(x, con_hifi.State_Solve(z_lofi), 'LineWidth', 3, 'color', 'black');
xlabel('$x$', 'Interpreter', 'latex');
ylabel('State');
legend({'$T$', '$\tilde{S}(\tilde{z})$', '$S(\tilde{z})$'}, 'Interpreter', 'latex', 'FontSize', 18, 'Location', 'south');
set(gca, 'fontsize', 18);
cd(write_path);
saveas(gca, 'states_at_lofi_source', 'epsc');
cd(working_path);

figure;
hold on;
plot(x, prior_delta_samples, 'LineWidth', 3, 'color', [.9, .9, .9]);
plot(x, prior_delta_samples(:, 1:10), 'LineWidth', 3);
xlabel('$x$', 'Interpreter', 'latex');
ylabel('Prior $\delta(\tilde{z},\theta)$', 'Interpreter', 'latex');
set(gca, 'fontsize', 18);
cd(write_path);
saveas(gca, 'prior_state_samples', 'epsc');
cd(working_path);

figure;
hold on;
plot(x, z_samples, 'LineWidth', 3);
xlabel('$x$', 'Interpreter', 'latex');
ylabel('');
yticks([]);
set(gca, 'fontsize', 18);
cd(write_path);
saveas(gca, 'prior_z_samples', 'epsc');
cd(working_path);

figure;
hold on;
plot(x, Z_test(:, 1), 'LineWidth', 3);
plot(x, Z_test(:, 2), 'LineWidth', 3);
plot(x, Z_test(:, 3), 'LineWidth', 3);
xlabel('$x$', 'Interpreter', 'latex');
ylabel('Source');
legend({'$z_1$', '$z_2$', '$z_{rep}$'}, 'Interpreter', 'latex');
set(gca, 'fontsize', 18);
cd(write_path);
saveas(gca, 'z_data', 'epsc');
cd(working_path);

figure;
hold on;
plot(x, md_update.post_data.D(:, 1), 'color', 'cyan', 'LineWidth', 3);
plot(x, delta_mean{1}, '--', 'color', 'red', 'LineWidth', 3);
for k = 1:num_post_samples
    plot(x, delta_samples{1}(:, k), 'color', [.9, .9, .9], 'LineWidth', 3);
end
plot(x, md_update.post_data.D(:, 1), 'color', 'cyan', 'LineWidth', 3);
plot(x, delta_mean{1}, '--', 'color', 'red', 'LineWidth', 3);
xlabel('$x$', 'Interpreter', 'latex');
ylabel('Discrepancy');
ylim([-3, 5]);
legend({'$d_1$', '$\delta(z_1,\theta)$'}, 'Interpreter', 'latex', 'Location', 'northwest');
set(gca, 'fontsize', 18);
cd(write_path);
saveas(gca, 'post_discrepancy_1', 'epsc');
cd(working_path);

figure;
hold on;
plot(x, md_update.post_data.D(:, 2), 'color', 'cyan', 'LineWidth', 3);
plot(x, delta_mean{2}, '--', 'color', 'red', 'LineWidth', 3);
for k = 1:num_post_samples
    plot(x, delta_samples{2}(:, k), 'color', [.9, .9, .9], 'LineWidth', 3);
end
plot(x, md_update.post_data.D(:, 2), 'color', 'cyan', 'LineWidth', 3);
plot(x, delta_mean{2}, '--', 'color', 'red', 'LineWidth', 3);
xlabel('$x$', 'Interpreter', 'latex');
ylabel('Discrepancy');
ylim([-3, 5]);
legend({'$d_2$', '$\delta(z_2,\theta)$'}, 'Interpreter', 'latex', 'Location', 'northwest');
set(gca, 'fontsize', 18);
cd(write_path);
saveas(gca, 'post_discrepancy_2', 'epsc');
cd(working_path);

figure;
hold on;
plot(x, delta_mean{3}, '--', 'color', 'red', 'LineWidth', 3);
for k = 1:num_post_samples
    plot(x, delta_samples{3}(:, k), 'color', [.9, .9, .9], 'LineWidth', 3);
end
plot(x, delta_mean{3}, '--', 'color', 'red', 'LineWidth', 3);
xlabel('$x$', 'Interpreter', 'latex');
ylabel('Discrepancy');
legend('$\delta(z_{ref},\theta)$', 'Interpreter', 'latex', 'Location', 'northwest');
set(gca, 'fontsize', 18);
cd(write_path);
saveas(gca, 'post_discrepancy_ref', 'epsc');
cd(working_path);

%%
L = length(z_update_mean_range);
mean_error = zeros(L, 1);
post_var = zeros(L, 1);
normalize = sqrt(z_hifi' * md_interface.Apply_M_z(z_hifi));
for k = 1:L
    e = z_hifi - z_update_mean_range{k};
    mean_error(k) = sqrt(e' * md_interface.Apply_M_z(e)) / normalize;

    C = cov(z_update_samples_range{k}');
    post_var(k) = ones(m, 1)' * md_interface.Apply_M_z(diag(C));
end

figure;
hold on;
yyaxis left;
plot(rank_range, mean_error, 'LineWidth', 3);
ylabel('Mean Update Relative Error');
yyaxis right;
plot(rank_range, post_var, 'LineWidth', 3);
ylabel('Optimal Solution Posterior Variance');
xlabel('Subspace Rank');
set(gca, 'fontsize', 18);
cd(write_path);
saveas(gca, 'vary_rank', 'epsc');
cd(working_path);

figure;
semilogy(md_update.gevp.evals / md_update.gevp.evals(1), 'LineWidth', 3, 'color', 'black');
xlabel('Subspace Rank');
ylabel('Normalized Generalized Eigenvalue');
set(gca, 'fontsize', 18);
cd(write_path);
saveas(gca, 'gen_eigs', 'epsc');
cd(working_path);

%%
j = 4;
u_true_update_4 = con_hifi.State_Solve(z_update_mean_range{j});
Jhat_update_4 = opt_hifi.Jhat(z_update_mean_range{j});

u_true_update_samples_4 = zeros(m, num_post_samples);
Jhat_update_samples_4 = zeros(num_post_samples, 1);
for k = 1:num_post_samples
    u_true_update_samples_4(:, k) = con_hifi.State_Solve(z_update_samples_range{j}(:, k));
    Jhat_update_samples_4(k) = opt_hifi.Jhat(z_update_samples_range{j}(:, k));
end

figure;
hold on;
plot(x, z_lofi, 'color', 'black', 'LineWidth', 3);
plot(x, z_hifi, 'color', 'cyan', 'LineWidth', 3);
plot(x, z_update_mean_range{j}, '--', 'color', 'red', 'LineWidth', 3);
for k = 1:num_post_samples
    plot(x, z_update_samples_range{j}(:, k), 'color', [.9, .9, .9], 'LineWidth', 3);
end
plot(x, md_update.z_opt, 'color', 'black', 'LineWidth', 3);
plot(x, z_hifi, 'color', 'cyan', 'LineWidth', 3);
plot(x, z_update_mean_range{j}, '--', 'color', 'red', 'LineWidth', 3);
legend({'$\tilde{z}$', '$z^\star$', '$\overline{z}$'}, 'Location', 'south', 'Interpreter', 'latex');
ylim([-100, 400]);
set(gca, 'fontsize', 18);
cd(write_path);
saveas(gca, 'post_opt_solution_rank_4', 'epsc');
cd(working_path);

figure;
hold on;
histogram(Jhat_update_samples_4);
plot([Jhat_update_4, Jhat_update_4], [0, 80], '--', 'LineWidth', 3, 'Color', 'red');
plot([Jhat_hifi, Jhat_hifi], [0, 80], 'LineWidth', 3, 'Color', 'cyan');
plot([Jhat_lofi, Jhat_lofi], [0, 80], 'LineWidth', 3, 'Color', 'black');
xlabel('High-fidelity objective function value');
yticks([]);
set(gca, 'fontsize', 18);
cd(write_path);
saveas(gca, 'post_obj_fun_vals_rank_4', 'epsc');
cd(working_path);

%%
j = 11;
u_true_update_11 = con_hifi.State_Solve(z_update_mean_range{j});
Jhat_update_11 = opt_hifi.Jhat(z_update_mean_range{j});

u_true_update_samples_11 = zeros(m, num_post_samples);
Jhat_update_samples_11 = zeros(num_post_samples, 1);
for k = 1:num_post_samples
    u_true_update_samples_11(:, k) = con_hifi.State_Solve(z_update_samples_range{j}(:, k));
    Jhat_update_samples_11(k) = opt_hifi.Jhat(z_update_samples_range{j}(:, k));
end

figure;
hold on;
plot(x, z_lofi, 'color', 'black', 'LineWidth', 3);
plot(x, z_hifi, 'color', 'cyan', 'LineWidth', 3);
plot(x, z_update_mean_range{j}, '--', 'color', 'red', 'LineWidth', 3);
for k = 1:num_post_samples
    plot(x, z_update_samples_range{j}(:, k), 'color', [.9, .9, .9], 'LineWidth', 3);
end
plot(x, md_update.z_opt, 'color', 'black', 'LineWidth', 3);
plot(x, z_hifi, 'color', 'cyan', 'LineWidth', 3);
plot(x, z_update_mean_range{j}, '--', 'color', 'red', 'LineWidth', 3);
legend({'$\tilde{z}$', '$z^\star$', '$\overline{z}$'}, 'Location', 'south', 'Interpreter', 'latex');
ylim([-100, 400]);
set(gca, 'fontsize', 18);
cd(write_path);
saveas(gca, 'post_opt_solution_rank_11', 'epsc');
cd(working_path);

figure;
hold on;
histogram(Jhat_update_samples_11);
plot([Jhat_update_11, Jhat_update_11], [0, 80], '--', 'LineWidth', 3, 'Color', 'red');
plot([Jhat_hifi, Jhat_hifi], [0, 80], 'LineWidth', 3, 'Color', 'cyan');
plot([Jhat_lofi, Jhat_lofi], [0, 80], 'LineWidth', 3, 'Color', 'black');
xlabel('High-fidelity objective function value');
yticks([]);
set(gca, 'fontsize', 18);
cd(write_path);
saveas(gca, 'post_obj_fun_vals_rank_11', 'epsc');
cd(working_path);
