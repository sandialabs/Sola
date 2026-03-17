%%
clear;
close all;
rng(121234);

suppress_figures = true;

m = 51;
x = linspace(0, 1, m)';
[M, S] = Assemble_Mass_and_Stiffness(m);

data_interface = MD_Data_Interface_synthetic_test_control_vec(m);
opt_prob_interface = MD_Opt_Prob_Interface_synthetic_test_control_vec(m);

u_hyperparam_interface = MD_u_Hyperparameter_Interface_synthetic_test_control_vec(m);
u_prior_interface = MD_Numeric_Laplacian_u_Prior_Interface(S, M, data_interface, u_hyperparam_interface);

num_state_solves = 100;
M_z = eye(2);
M_z(1, 1) = 2;
z_hyperparam_interface = MD_z_Hyperparameter_Interface_synthetic_test_control_vec(num_state_solves, m);
z_prior_interface = MD_Vector_z_Prior_Interface(M_z, data_interface, z_hyperparam_interface, u_prior_interface);

%%
num_prior_samples = 100;
md_prior_sampling = MD_Prior_Sampling(data_interface, u_prior_interface, z_prior_interface);
md_prior_sampling.Generate_Prior_Discrepancy_Sample_Data(num_prior_samples);
if ~suppress_figures
    md_prior_vis = MD_Prior_Visualization(md_prior_sampling);
    md_prior_vis.Visualization_for_Prior_Discrepancy_at_z_opt(1);
    md_prior_vis.Visualization_for_Prior_Discrepancy_at_z_pert(1);
end

%%
md_post_sampling = MD_Posterior_Sampling(data_interface, u_prior_interface, z_prior_interface);
alpha_d = u_prior_interface.u_hyperparam_interface.alpha_d;
num_post_samples = 100;
md_post_sampling.Compute_Posterior_Data(alpha_d, num_post_samples);

%%
md_hessian_analysis = MD_Hessian_Analysis(opt_prob_interface, z_prior_interface);
md_hessian_analysis.Compute_Hessian_GEVP(data_interface.z_opt, 2, 0);
md_update = MD_Update(md_post_sampling, md_hessian_analysis);

[z_update_mean, z_update_samples] = md_update.Posterior_Update_Samples();

if ~suppress_figures
    z_hifi = [1 - 2 * data_interface.epsilon; 1 - data_interface.epsilon];
    figure;
    for k = 1:2
        subplot(2, 1, k);
        hold on;
        plot(z_update_samples(k, :), ones(num_post_samples, 1), 'o', 'Color', [.9, .9, .9], 'MarkerSize', 10);
        plot(z_update_mean(k), 1, 'o', 'color', 'blue', 'MarkerSize', 10);
        plot(z_hifi(k), 1, '*', 'color', 'red', 'MarkerSize', 10);
    end
end

%%
z_mean_ref = load('reference_solution.mat').z_update_mean;
z_samples_ref = load('reference_solution.mat').z_update_samples;
ref_diff = max(norm(z_mean_ref - z_update_mean) / norm(z_update_mean), norm(z_update_samples - z_samples_ref) / norm(z_update_samples));

if ref_diff > 1.e-9
    fprintf(2, '\nmodel_discrepancy/synthetic_test_control_vec failed.\n');
else
    fprintf(1, '\nmodel_discrepancy/synthetic_test_control_vec passed.\n');
end
