clear;
close all;
addpath(genpath('../../src/'));
clc;

n_y = 50;
n_z = n_y;

ballistic_data = load('~/Documents/alpha/examples/Ballistic_2DOF/Forward_Solve_Data.mat');

p = length(ballistic_data.f);
posterior_uncertainty = zeros(p, 1);

for k = 1:p
    T = ballistic_data.t{k}(end);
    n_t = length(ballistic_data.t{k});

    con = Thermal_Constraint(n_y, n_z, T, n_t);
    z_true = (1.e-2) * (2 + cos(3 * pi * con.x));

    forcing = @(x, t) interp1(ballistic_data.t{k}, ballistic_data.f{k}, t) * (x.^2);
    Generate_Obs_Data(con, z_true, forcing);

    con.forcing = load('Obs_Data.mat', 'forcing').forcing;
    likelihood = Thermal_Likelihood_Model(n_y, n_t);
    prior = Thermal_Prior_Model(con);
    bayes_inv = Bayesian_Inversion(likelihood, prior, con);

    z0 = load('Obs_Data.mat', 'z_true').z_true;
    bayes_inv.opt.z_lb = zeros(n_y, 1);
    bayes_inv.opt.max_cg_iter = 500;
    bayes_inv.opt.step_tol = 1.e-8;
    [u_map, z_map] = bayes_inv.Compute_MAP_Point(z0);

    [~, ~, hessian_data] = bayes_inv.opt.Jhat(z_map);
    H = bayes_inv.opt.Jhat_hessVec(hessian_data, eye(n_y));
    posterior_uncertainty(k) = trace(linsolve(H, con.M));
end
