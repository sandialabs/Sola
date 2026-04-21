%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
clear;
close all;
% rng(121234);

suppress_figures = false;

m = 51;

data_interface = MD_OUU_Data_Interface_synthetic_test_OUU();
data_interface.Load_Data();

Xi = load('Optimization_Results.mat', 'Xi').Xi;
N = size(Xi, 2);
obj = Synthetic_Test_OUU_Objective(m);
cons = cell(N, 1);
for k = 1:N
    cons{k} = Synthetic_Test_OUU_Constraint(Xi(:, k));
end
opt = Reduced_Space_Optimization_Under_Uncertainty(obj, cons);
x = obj.x;

opt_prob_interface = MD_OUU_Opt_Prob_Interface_Sola(data_interface, opt);
us_prior_interface = MD_u_Prior_Interface_synthetic_test_OUU(m);

ensemble_weighting = MD_OUU_Ensemble_Weighting_Matrix(data_interface, us_prior_interface);
u_prior_interface = MD_OUU_u_Prior_Interface(us_prior_interface, data_interface, ensemble_weighting);
z_prior_interface = MD_z_Prior_Interface_synthetic_test_OUU(m);

%%
num_prior_samples = 100;
md_prior_sampling = MD_Prior_Sampling(data_interface, u_prior_interface, z_prior_interface);
prior_delta_samples_zopt = md_prior_sampling.Prior_Discrepancy_Samples_at_z_opt(num_prior_samples);

if ~suppress_figures
    I = 1:N:(m * N);
    u = prior_delta_samples_zopt(I, :);

    figure;
    plot(x, u, 'LineWidth', 3);
end

if ~suppress_figures
    k = 1;
    u = data_interface.Reshape_State_to_Mat(prior_delta_samples_zopt(:, k));

    c = ensemble_weighting.C(1, :);
    c_normalized = (c - min(c)) / (max(c) - min(c));
    n = size(u, 2);
    figure;
    hold on;
    col = colormap;
    for j = 1:n
        color = col(round(c_normalized(j) * (size(colormap, 1) - 1)) + 1, :);
        plot(x, u(:, j), 'Color', color, 'LineWidth', 3);
    end
    hold off;
    colorbar;
    caxis([min(c) max(c)]);
end

%%
x = opt.obj.x;
z = zeros(m, 3);
z(:, 1) = x;
z(:, 2) = x.^2 + 1;
z(:, 3) = sin(2 * pi * x);
prior_delta_samples_z = md_prior_sampling.Prior_Discrepancy_Samples(z, num_prior_samples);

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
    i = 2;

    k = 20;
    u = data_interface.Reshape_State_to_Mat(delta_samples{i}(:, k));
    c = ensemble_weighting.C(1, :);
    c_normalized = (c - min(c)) / (max(c) - min(c));
    n = size(u, 2);
    figure;
    hold on;
    col = colormap;
    for j = 1:n
        color = col(round(c_normalized(j) * (size(colormap, 1) - 1)) + 1, :);
        plot(x, u(:, j), 'Color', color, 'LineWidth', 3);
    end
    hold off;
    colorbar;
    caxis([min(c) max(c)]);

    u1 = data_interface.Reshape_State_to_Mat(data_interface.D(:, i));
    u2 = data_interface.Reshape_State_to_Mat(delta_mean{i});
    n = size(u1, 2);
    figure;
    hold on;
    col = lines(n);
    for j = 1:n
        color = col(j, :);
        plot(x, u1(:, j), 'Color', color, 'LineWidth', 3);
        plot(x, u2(:, j), '--', 'Color', 'black', 'LineWidth', 3);
    end
    hold off;
end

if ~suppress_figures
    i = 3;

    k = 20;
    u = data_interface.Reshape_State_to_Mat(delta_samples{i}(:, k));
    c = ensemble_weighting.C(1, :);
    c_normalized = (c - min(c)) / (max(c) - min(c));
    n = size(u, 2);
    figure;
    hold on;
    col = colormap;
    for j = 1:n
        color = col(round(c_normalized(j) * (size(colormap, 1) - 1)) + 1, :);
        plot(x, u(:, j), 'Color', color, 'LineWidth', 3);
    end
    hold off;
    colorbar;
    caxis([min(c) max(c)]);
end

%%
md_hessian_analysis = MD_Hessian_Analysis(opt_prob_interface, z_prior_interface);
num_evals = 8;
oversampling = 10;
md_hessian_analysis.Compute_Hessian_GEVP(data_interface.z_opt, num_evals, oversampling);

md_update = MD_Update(md_post_sampling, md_hessian_analysis);
[z_update_mean, z_update_samples] = md_update.Posterior_Update_Samples();

%%
if ~suppress_figures
    z_hifi = load('Optimization_Results.mat', 'z_hifi').z_hifi;
    figure;
    hold on;
    plot(x, data_interface.z_opt, 'Color', 'magenta', 'LineWidth', 3);
    plot(x, z_hifi, 'Color', 'cyan', 'LineWidth', 3);
    plot(x, z_update_mean, '--', 'Color', 'red', 'LineWidth', 3);
end
