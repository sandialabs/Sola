clear;
close all;
addpath(genpath('../../../src'));
rng(1342);

% Instantiate the Example_1 object
con = Example_1_Constraint();
likelihood = Likelihood_Model_Example_1();
prior = Prior_Model_Example_1();
bayes_inv = Bayesian_Inversion(likelihood, prior, con);

z0 = rand(2, 1);
bayes_inv.opt.verbose = false;
[u_map, z_map] = bayes_inv.Compute_MAP_Point(z0);

u_sol = load('Solution_Example_1.mat', 'u_map').u_map;
z_sol = load('Solution_Example_1.mat', 'z_map').z_map;

error = 0;
error = max(error, norm(u_sol - u_map));
error = max(error, norm(z_sol - z_map));

if error > 1.e-11
    fprintf(2,'\nbayesian_inversion/Example_1 failed.\n');
else
    fprintf(1,'\nbayesian_inversion/Example_1 passed.\n');
end

% save('Solution_Example_1.mat','u_map','z_map')
