%%
clear;
close all;
addpath(genpath('../../../src'));
%rng(121234);

suppress_figures = false;

m = 51;

data_interface = MD_OUU_Data_Interface_synthetic_test_OUU();
data_interface.Load_Data();

Xi = data_interface.Xi;
N = size(Xi, 2);
obj = Synthetic_Test_OUU_Objective(m);
cons = cell(N, 1);
for k = 1:N
    cons{k} = Synthetic_Test_OUU_Constraint(Xi(:, k));
end
opt = Reduced_Space_Optimization_Under_Uncertainty(obj, cons);
x = obj.x;

opt_prob_interface = MD_OUU_Opt_Prob_Interface_Sabl(data_interface, opt);

us_prior_interface = MD_u_Prior_Interface_synthetic_test_OUU(m);
u_prior_interface = MD_OUU_u_Prior_Interface(us_prior_interface, data_interface);
z_prior_interface = MD_z_Prior_Interface_synthetic_test_OUU(m);

%%
num_prior_samples = 100;
md_prior_sampling = MD_Prior_Sampling(data_interface, u_prior_interface, z_prior_interface);
prior_delta_samples_zopt = md_prior_sampling.Prior_Discrepancy_Samples_at_z_opt(num_prior_samples);

if ~suppress_figures
    k = 1;
    u = data_interface.Reshape_State_to_Mat(prior_delta_samples_zopt(:,k));

    colormap('jet');
    c = u_prior_interface.L(k,:);
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
    clim([min(c) max(c)]); 
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

%%
md_hessian_analysis = MD_Hessian_Analysis(opt_prob_interface, z_prior_interface);
num_evals = 20;
oversampling = 10;
md_hessian_analysis.Compute_Hessian_GEVP(data_interface.z_opt, num_evals, oversampling);

md_update = MD_Update(md_post_sampling, md_hessian_analysis);
[z_update_mean, z_update_samples] = md_update.Posterior_Update_Samples();
