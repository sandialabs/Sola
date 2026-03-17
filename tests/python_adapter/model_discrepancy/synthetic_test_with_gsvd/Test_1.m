%%
clear;
close all;
rng(121234);

update_python = false;

if update_python
    clear classes;
    update_python = true;
end

suppress_figures = true;

m = int64(51);
x = linspace(0, 1, m)';

alpha_u = 1 / (2^2);
alpha_z = 1 / (100^2);

mod = py.importlib.import_module('MD_Data_Interface_Python_Synthetic_Test');
if update_python
    py.importlib.reload(mod);
end
data_interface_python = mod.MD_Data_Interface_Python_Synthetic_Test(m);
data_interface = MD_Data_Interface_Py(data_interface_python);
data_interface.Load_Data();

mod = py.importlib.import_module('MD_Opt_Prob_Interface_Python_Synthetic_Test');
if update_python
    py.importlib.reload(mod);
end
opt_prob_interface_python = mod.MD_Opt_Prob_Interface_Python_Synthetic_Test(m);
opt_prob_interface = MD_Opt_Prob_Interface_Py(opt_prob_interface_python);

mod = py.importlib.import_module('MD_Elliptic_u_Prior_Interface_Python_Synthetic_Test');
if update_python
    py.importlib.reload(mod);
end
u_prior_interface_python = mod.MD_Elliptic_u_Prior_Interface_Python_Synthetic_Test(m);
u_prior_interface = MD_Elliptic_u_Prior_Interface_Py(u_prior_interface_python, alpha_u);

mod = py.importlib.import_module('MD_Elliptic_z_Prior_Interface_Python_Synthetic_Test');
if update_python
    py.importlib.reload(mod);
end
z_prior_interface_python = mod.MD_Elliptic_z_Prior_Interface_Python_Synthetic_Test(m, alpha_z);
z_prior_interface = MD_Elliptic_z_Prior_Interface_Py(z_prior_interface_python, alpha_z);

num_sing_vals = 50;
oversampling = 1;
num_subspace_iters = 2;
u_vec = zeros(m, 1);
u_prior_interface.Compute_E_u_Inverse_GSVD(num_sing_vals, oversampling, num_subspace_iters, u_vec);

num_prior_samples = 100;
md_prior_sampling = MD_Prior_Sampling(data_interface, u_prior_interface, z_prior_interface);

%%
delta_samples = md_prior_sampling.Prior_Discrepancy_Samples_at_z_opt(num_prior_samples);

if ~suppress_figures
    figure;
    plot(x, delta_samples, 'LineWidth', 3, 'color', [.9, .9, .9]);
end

z = zeros(m, 3);
z(:, 1) = x;
z(:, 2) = x.^2 + 1;
z(:, 3) = sin(2 * pi * x);
delta_prior_samples = md_prior_sampling.Prior_Discrepancy_Samples(z, num_prior_samples);
if ~suppress_figures
    for k = 1:10
        figure;
        hold on;
        plot(x, delta_prior_samples{k}, 'LineWidth', 3);
    end
end

%%
md_post_sampling = MD_Posterior_Sampling(data_interface, u_prior_interface, z_prior_interface);
alpha_d = 1.e-5;
num_post_samples = 100;
md_post_sampling.Compute_Posterior_Data(alpha_d, num_post_samples);
Z_test = randn(m, 3);
Z_test(:, 1:2) = md_post_sampling.post_data.Z;
Z_test(:, 3) = 1.5 * ones(m, 1);
[delta_mean, delta_samples] = md_post_sampling.Posterior_Discrepancy_Samples(Z_test);

if ~suppress_figures
    figure;
    hold on;
    plot(x, md_post_sampling.post_data.D(:, 1), 'color', 'black', 'LineWidth', 3);
    plot(x, delta_mean{1}, '--', 'color', 'red', 'LineWidth', 3);
    for k = 1:num_post_samples
        plot(x, delta_samples{1}(:, k), 'color', [.9, .9, .9], 'LineWidth', 3);
    end
    plot(x, md_post_sampling.post_data.D(:, 1), 'color', 'black', 'LineWidth', 3);
    plot(x, delta_mean{1}, '--', 'color', 'red', 'LineWidth', 3);

    figure;
    hold on;
    plot(x, md_post_sampling.post_data.D(:, 2), 'color', 'black', 'LineWidth', 3);
    plot(x, delta_mean{2}, '--', 'color', 'red', 'LineWidth', 3);
    for k = 1:num_post_samples
        plot(x, delta_samples{2}(:, k), 'color', [.9, .9, .9], 'LineWidth', 3);
    end
    plot(x, md_post_sampling.post_data.D(:, 2), 'color', 'black', 'LineWidth', 3);
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
md_hessian_analysis = MD_Hessian_Analysis(opt_prob_interface, z_prior_interface);
md_update = MD_Update(md_post_sampling, md_hessian_analysis);
[z_update_mean, z_update_samples] = md_update.Posterior_Update_Samples();

if ~suppress_figures
    figure;
    hold on;
    plot(x, (1 + x) / (1.2^(1 / 3)), 'color', 'black', 'LineWidth', 3);
    plot(x, 1 + x, 'color', 'cyan', 'LineWidth', 3);
    plot(x, z_update_mean, '--', 'color', 'red', 'LineWidth', 3);
    for k = 1:num_post_samples
        plot(x, z_update_samples(:, k), 'color', [.9, .9, .9], 'LineWidth', 3);
    end
    plot(x, (1 + x) / (1.2^(1 / 3)), 'color', 'black', 'LineWidth', 3);
    plot(x, 1 + x, 'color', 'cyan', 'LineWidth', 3);
    plot(x, z_update_mean, '--', 'color', 'red', 'LineWidth', 3);
end

%%
% save('reference_solution.mat','z_update_mean')

z_mean_ref = load('reference_solution.mat').z_update_mean;
ref_diff = norm(z_mean_ref - z_update_mean);

if ref_diff > 1.e-9
    fprintf(2, '\npython_adapter/model_discrepancy/synthetic_test_with_gsvd failed.\n');
else
    fprintf(1, '\npython_adapter/model_discrepancy/synthetic_test_with_gsvd passed.\n');
end
