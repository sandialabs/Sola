clear;
close all;
clc;
addpath(genpath('../../src'));

m = 50;
con = Thermal_Constraint(m);
con.forcing = load('Obs_Data.mat', 'forcing').forcing;
likelihood = Thermal_Likelihood_Model(m);
prior = Thermal_Prior_Model(con);
bayes_inv = Bayesian_Inversion(likelihood, prior, con);

plot_prior_samples = true;
if plot_prior_samples
    num_samps = 20;
    Z_prior = prior.Compute_Prior_Samples(num_samps);
    plot(con.x, Z_prior);
end

z0 = prior.Get_Prior_Mean();
bayes_inv.opt.z_lb = zeros(m, 1);
bayes_inv.opt.max_cg_iter = 500;
[u_map, z_map] = bayes_inv.Compute_MAP_Point(z0);

figure;
hold on;
u_true = load('Obs_Data.mat', 'u_true').u_true;
u_data = load('Obs_Data.mat', 'u_data').u_data;
plot(con.x, u_map, 'LineWidth', 3);
plot(con.x, u_true, '--', 'LineWidth', 3);
plot(con.x(likelihood.obs_vec), u_data(likelihood.obs_vec), 'o', 'MarkerSize', 10);
title('State');
legend({'MAP Estimate', 'Truth'});
set(gca, 'fontsize', 18);

figure;
hold on;
z_true = load('Obs_Data.mat', 'z_true').z_true;
plot(con.x, z_map, 'LineWidth', 3);
plot(con.x, z_true, '--', 'LineWidth', 3);
title('Conductivity');
legend({'MAP Estimate', 'Truth'});
set(gca, 'fontsize', 18);
