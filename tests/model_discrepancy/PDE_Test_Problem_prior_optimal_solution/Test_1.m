%%
clear;
close all;
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

z_tilde = load('z_opt.mat').z_opt;
z_star = load('z_hifi.mat').z_hifi;

%%
data_interface = MD_Data_Interface_PDE_Test_Problem();

alpha_u = 1 / (2^2);
alpha_z = 1 / (3^2);
u_prior_interface = MD_Elliptic_u_Prior_Interface_PDE_Test_Problem(alpha_u, opt_lofi);
z_prior_interface = MD_Elliptic_z_Prior_Interface_PDE_Test_Problem(alpha_z, opt_lofi);

opt_prob_interface = MD_Opt_Prob_Interface_Sabl(opt_lofi, data_interface);
md_hessian_analysis = MD_Hessian_Analysis(opt_prob_interface, z_prior_interface);

num_evals = 10;
oversampling = 10;
md_hessian_analysis.Compute_Hessian_GEVP(z_tilde, num_evals, oversampling);

%%
num_prior_samples = 100;
md_prior_opt_sol_sampling = MD_Prior_Opt_Solution_Sampling(data_interface, u_prior_interface, z_prior_interface, opt_prob_interface, md_hessian_analysis);

z_star_prior_samps = md_prior_opt_sol_sampling.Generate_Prior_Opt_Solution_Samples(num_prior_samples);

if ~suppress_figures
    figure;
    hold on;
    plot(x, z_tilde, 'color', 'black', 'LineWidth', 3);
    plot(x, z_star, '--', 'color', 'red', 'LineWidth', 3);
    plot(x, z_star_prior_samps, 'Color', .9 * ones(3, 1), 'LineWidth', 3);
    plot(x, z_tilde, 'color', 'black', 'LineWidth', 3);
    plot(x, z_star, '--', 'color', 'red', 'LineWidth', 3);
    legend({'Low-fidelity', 'High-fidelity'});
    set(gca, 'fontsize', 18);
end

%%
ref = load('reference_solution.mat');
ref_diff = norm(z_star_prior_samps - ref.z_star_prior_samps, 'fro');

if ref_diff > 1.e-14
    fprintf(2, '\nmodel_discrepancy/PDE_Test_Problem_prior_optimal_solution failed.\n');
else
    fprintf(1, '\nmodel_discrepancy/PDE_Test_Problem_prior_optimal_solution passed.\n');
end
