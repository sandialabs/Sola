clear;
close all;
clc;
addpath(genpath('../../src/'));
rng(154);

write_figure_to_file = false;
working_path = pwd;
write_path = '~/Desktop';

m = 100;
theta = 1 + (1 - linspace(0, 1, m)').^2;
con = Adv_Diff_Constraint(theta);
Generate_Obs_Data(con);
likelihood = Adv_Diff_Likelihood_Model(m);
prior = Adv_Diff_Prior_Model(con);
bayes_inv = Bayesian_Inversion(likelihood, prior, con);
bayes_inv.opt.fun_tol = 1.e-10;
bayes_inv.opt.opt_tol = 1.e-10;

plot_prior = false;
if plot_prior
    num_samps = 20;
    z_prior_samples = prior.Compute_Prior_Samples(num_samps);
    figure;
    plot(con.x, z_prior_samples, 'LineWidth', 3);
    set(gca, 'fontsize', 18);
end

z0 = prior.Get_Prior_Mean();
verbose = 'iter';
[u_opt, z_opt] = bayes_inv.Compute_MAP_Point(z0);

z_true = load('Obs_Data.mat', 'z_true').z_true;
figure;
hold on;
plot(con.x, z_true, 'LineWidth', 3);
plot(con.x, z_opt, '--', 'LineWidth', 3);
legend({'$z^\dagger$', '$\overline{z}$'}, 'Location', 'best', 'Interpreter', 'latex');
xlabel('$x$', 'Interpreter', 'latex');
ylabel('Diffusion Coefficient');
set(gca, 'fontsize', 20);
if write_figure_to_file
    cd(write_path);
    saveas(gca, 'map_point', 'epsc');
    cd(working_path);
end

u_true = load('Obs_Data.mat', 'u_true').u_true;
obs_points = bayes_inv.likelihood.Observation_Operator_Apply(con.x);
u_data = likelihood.Get_Observed_Data();
figure;
hold on;
plot(con.x, u_true, 'LineWidth', 3);
plot(con.x, u_opt, '--', 'LineWidth', 3);
scatter(obs_points, u_data, 100, 'filled');
legend({'$u(z^\dagger,\theta^\dagger)$', '$u(\overline{z},\overline{\theta})$', 'Observed Data'}, 'Location', 'best', 'Interpreter', 'latex');
xlabel('$x$', 'Interpreter', 'latex');
ylabel('State');
set(gca, 'fontsize', 20);
if write_figure_to_file
    cd(write_path);
    saveas(gca, 'state_solution', 'epsc');
    cd(working_path);
end

save('Optimization_Results.mat', 'z_opt', 'theta');
