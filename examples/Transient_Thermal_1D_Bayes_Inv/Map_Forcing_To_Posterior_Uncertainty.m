clear;
close all;
clc;

m = 50;
n = m;
T = 10;
N = 100;

con = Thermal_Constraint(m, n, T, N);
z_true = (1.e-2) * (2 + cos(3 * pi * con.x));

forcings = cell(4, 1);
forcings{1} = @(x, t) (1 + 0 * t) * exp(-100 * (x - 0.5).^2);
forcings{2} = @(x, t) (1 + .5 * sin(2 * pi * t / T)) * exp(-100 * (x - 0.5).^2);
forcings{3} = @(x, t) (1 + t / T) * exp(-100 * (x - 0.5).^2);
forcings{4} = @(x, t) (1 + t.^2 / T^2) * exp(-100 * (x - 0.5).^2);

p = length(forcings);
posterior_uncertainty = zeros(p, 1);

for k = 1:p
    forcing = forcings{k};
    Generate_Obs_Data(con, z_true, forcing);

    con.forcing = load('Obs_Data.mat', 'forcing').forcing;
    likelihood = Thermal_Likelihood_Model(m, N);
    prior = Thermal_Prior_Model(con);
    bayes_inv = Bayesian_Inversion(likelihood, prior, con);

    z0 = load('Obs_Data.mat', 'z_true').z_true;
    bayes_inv.opt.z_lb = zeros(m, 1);
    bayes_inv.opt.max_cg_iter = 500;
    bayes_inv.opt.step_tol = 1.e-8;
    [u_map, z_map] = bayes_inv.Compute_MAP_Point(z0);

    [~, ~, hessian_data] = bayes_inv.opt.Jhat(z_map);
    H = bayes_inv.opt.Jhat_hessVec(hessian_data, eye(m));
    posterior_uncertainty(k) = trace(linsolve(H, con.M));
end
