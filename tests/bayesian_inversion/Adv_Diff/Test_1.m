clear;
close all;
rng(1423435);

%%%%% Create on instance of this problem under Bayesian inversion and another under OED %%%%%

m = 200;
diff_coeff = 2;
vel_coeff = 10;
robin_coeff = 3;
con = Adv_Diff_Constraint(m, diff_coeff, vel_coeff, robin_coeff);

likelihood = Adv_Diff_Likelihood_Model(m);
prior = Adv_Diff_Prior_Model(con);

num_samps = 20;
Z_prior = prior.Compute_Prior_Samples(num_samps);

bayes_inv = Bayesian_Inversion(likelihood, prior, con);
bayes_inv.opt.verbose = false;

z0 = rand(m, 1);
[u_map, z_map] = bayes_inv.Compute_MAP_Point(z0);

u_sol = load('Solution_Adv_Diff.mat', 'u_map').u_map;
z_sol = load('Solution_Adv_Diff.mat', 'z_map').z_map;
Z_sol = load('Solution_Adv_Diff.mat', 'Z_prior').Z_prior;

error = 0;
error = max(error, norm(u_sol - u_map));
error = max(error, norm(z_sol - z_map));
error = max(error, norm(Z_sol - Z_prior));
if error > 1.e-11
    fprintf(2,'\nbayesian_inversion/Adv_Diff failed.\n');
else
    fprintf(1,'\nbayesian_inversion/Adv_Diff passed.\n');
end
% save('Solution_Adv_Diff.mat','u_map','z_map','Z_prior')
