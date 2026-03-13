clear;
close all;
addpath(genpath('../../src'));
rng(1423435);

m = 50;
diff_coeff = 1;
con = Poisson_Constraint(m, diff_coeff);
x = con.x;

likelihood = Poisson_Likelihood_Model(con);
prior = Poisson_Prior_Model(con);
bayes_inv = Bayesian_Inversion(likelihood, prior, con);

z0 = rand(m, 1);
bayes_inv.opt.verbose = false;
[u_map, z_map] = bayes_inv.Compute_MAP_Point(z0);

u_sol = load('Solution_Poisson.mat', 'u_map').u_map;
z_sol = load('Solution_Poisson.mat', 'z_map').z_map;

error = 0;
error = max(error, norm(u_sol - u_map));
error = max(error, norm(z_sol - z_map));
if error > 1.e-11
    fprintf(2,'\nbayesian_inversion/Poisson failed.\n');
else
    fprintf(1,'\nbayesian_inversion/Poisson  passed.\n');
end


% save('Solution_Poisson.mat','u_map','z_map')
