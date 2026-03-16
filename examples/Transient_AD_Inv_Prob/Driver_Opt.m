clear;
close all;
clc;

n_y = 50;
n_z = n_y;
T = 1;
n_t = 100;

con = Transient_Inv_Prob_Constraint_AD(n_y, n_z, T, n_t);
con.forcing = load('Obs_Data.mat', 'forcing').forcing;
con.AD_Initialization();
likelihood = Transient_Inv_Prob_Likelihood_Model(n_y, n_t);
prior = Transient_Inv_Prob_Prior_Model(con);
bayes_inv = Bayesian_Inversion(likelihood, prior, con);

plot_prior_samples = false;
if plot_prior_samples
    num_samps = 20;
    Z_prior = prior.Compute_Prior_Samples(num_samps);
    plot(con.x, Z_prior);
end

% z0 = prior.Get_Prior_Mean();
z0 = load('Obs_Data.mat', 'z_true').z_true;
bayes_inv.opt.z_lb = zeros(n_z, 1);
bayes_inv.opt.max_cg_iter = 500;
[u_map, z_map] = bayes_inv.Compute_MAP_Point(z0);

figure;
hold on;
z_true = load('Obs_Data.mat', 'z_true').z_true;
plot(con.x, z_map, 'LineWidth', 3);
plot(con.x, z_true, '--', 'LineWidth', 3);
title('Conductivity');
legend({'MAP Estimate', 'Truth'});
set(gca, 'fontsize', 18);
