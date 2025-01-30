%%
clear;
close all;
addpath(genpath('../../../src'));
rng(1234423);

suppress_figures = true;

m = 200;
diff_coeff = 1;
vel_coeff = 1 / 2;
robin_coeff = 2;
reg_coeff = 10;
obj = Adv_Diff_Objective(m, reg_coeff);
con_hifi = Adv_Diff_Constraint(m, diff_coeff, vel_coeff, robin_coeff);
con_lofi = Diff_Constraint(obj, con_hifi);
opt_hifi = Reduced_Space_Optimization(obj, con_hifi);
opt_lofi = Reduced_Space_Optimization(obj, con_lofi);
x = con_hifi.x;

%%
data_interface = MD_Data_Interface_PDE_Test_Problem();
data_interface.Load_Data();

alpha_u = 1 / (2^2);
alpha_z = 1 / (3^2);
u_prior_interface = MD_Elliptic_u_Prior_Interface_PDE_Test_Problem(alpha_u, opt_lofi);
z_prior_interface = MD_Elliptic_z_Prior_Interface_PDE_Test_Problem(alpha_z, opt_lofi);

%%
num_prior_samples = 100;
md_prior_sampling = MD_Prior_Sampling(data_interface, u_prior_interface, z_prior_interface);

delta_samples = md_prior_sampling.Prior_Discrepancy_Samples_at_z_opt(num_prior_samples);
if ~suppress_figures
    figure;
    plot(x, delta_samples(:, 1:10), 'LineWidth', 3);

    figure;
    plot(x, delta_samples, 'LineWidth', 3, 'color', [.9, .9, .9]);
end

%%
z = zeros(m, 3);
z(:, 1) = 1 + sin(2 * pi * x);
z(:, 2) = 1 + cos(2 * pi * x);
z(:, 3) = 1 + sin(20 * pi * x);
if ~suppress_figures
    figure;
    plot(x, z, 'LineWidth', 3);
end

delta_prior_samples = md_prior_sampling.Prior_Discrepancy_Samples(z, num_prior_samples);
if ~suppress_figures
    for k = 1:10
        figure;
        hold on;
        plot(x, delta_prior_samples{k}, 'LineWidth', 3);
    end
end

%%
md_post_samples = MD_Posterior_Sampling(data_interface, u_prior_interface, z_prior_interface);
alpha_d = 1.e-5;
num_post_samples = 100;
md_post_samples.Compute_Posterior_Data(alpha_d, num_post_samples);
Z_test = randn(m, 3);
Z_test(:, 1:2) = md_post_samples.post_data.Z;
Z_test(:, 3) = 1.5 * ones(m, 1);
[delta_mean, delta_samples] = md_post_samples.Posterior_Discrepancy_Samples(Z_test);

if ~suppress_figures
    figure;
    hold on;
    plot(x, md_post_samples.post_data.D(:, 1), 'color', 'black', 'LineWidth', 3);
    plot(x, delta_mean{1}, '--', 'color', 'red', 'LineWidth', 3);
    for k = 1:num_post_samples
        plot(x, delta_samples{1}(:, k), 'color', [.9, .9, .9], 'LineWidth', 3);
    end
    plot(x, md_post_samples.post_data.D(:, 1), 'color', 'black', 'LineWidth', 3);
    plot(x, delta_mean{1}, '--', 'color', 'red', 'LineWidth', 3);

    figure;
    hold on;
    plot(x, md_post_samples.post_data.D(:, 2), 'color', 'black', 'LineWidth', 3);
    plot(x, delta_mean{2}, '--', 'color', 'red', 'LineWidth', 3);
    for k = 1:num_post_samples
        plot(x, delta_samples{2}(:, k), 'color', [.9, .9, .9], 'LineWidth', 3);
    end
    plot(x, md_post_samples.post_data.D(:, 2), 'color', 'black', 'LineWidth', 3);
    plot(x, delta_mean{2}, '--', 'color', 'red', 'LineWidth', 3);

    figure;
    hold on;
    plot(x, delta_mean{3}, '--', 'color', 'red', 'LineWidth', 3);
    for k = 1:num_post_samples
        plot(x, delta_samples{3}(:, k), 'color', [.9, .9, .9], 'LineWidth', 3);
    end
    plot(x, delta_mean{3}, '--', 'color', 'red', 'LineWidth', 3);

end

%%
opt_prob_interface = MD_Opt_Prob_Interface_Sabl(opt_lofi, data_interface);
md_hessian_analysis = MD_Hessian_Analysis(opt_prob_interface, z_prior_interface);

num_evals = 10;
oversampling = 10;
md_hessian_analysis.Compute_Hessian_GEVP(data_interface.z_opt, num_evals, oversampling);

md_update = MD_Update(md_post_samples, md_hessian_analysis);

[z_update_mean, z_update_samples] = md_update.Posterior_Update_Samples();

if ~suppress_figures
    z_hifi = load('z_hifi.mat').z_hifi;
    figure;
    hold on;
    plot(x, md_update.z_opt, 'color', 'black', 'LineWidth', 3);
    plot(x, z_hifi, 'color', 'cyan', 'LineWidth', 3);
    plot(x, z_update_mean, '--', 'color', 'red', 'LineWidth', 3);
    for k = 1:num_post_samples
        plot(x, z_update_samples(:, k), 'color', [.9, .9, .9], 'LineWidth', 3);
    end
    plot(x, md_update.z_opt, 'color', 'black', 'LineWidth', 3);
    plot(x, z_hifi, 'color', 'cyan', 'LineWidth', 3);
    plot(x, z_update_mean, '--', 'color', 'red', 'LineWidth', 3);
end

%%
z_mean_ref = load('reference_solution.mat').z_update_mean;
z_samples_ref = load('reference_solution.mat').z_update_samples;
ref_diff = max(norm(z_mean_ref - z_update_mean) / norm(z_update_mean), norm(z_update_samples - z_samples_ref) / norm(z_update_samples));
if ref_diff > 1.e-14
    disp('PDE_Test_Problem difference:');
    disp(ref_diff);
end
