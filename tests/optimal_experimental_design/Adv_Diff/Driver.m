clear;
addpath(genpath('../../../src'));
rng(1423435);

m = 200;
diff_coeff = 2;
vel_coeff = 10;
robin_coeff = 3;
con = Adv_Diff_Constraint(m, diff_coeff, vel_coeff, robin_coeff);

likelihood = Adv_Diff_Likelihood_Model(m);
prior = Adv_Diff_Prior_Model(con);

num_trace_samples = 1000;
reguarlization_coeff = 1.e-1;
linear_oed = Linear_OED(likelihood, prior, con, num_trace_samples, reguarlization_coeff);
linear_oed.verbose = false;

num_sing_vals = 50;
oversampling = 10;
num_subspace_iters = 1;
sing_vals = linear_oed.Compute_Forward_Operator_GSVD(num_sing_vals, oversampling, num_subspace_iters);

w = linear_oed.Optimize_Design();

w_sol = load('Solution_Adv_Diff.mat', 'w').w;
sing_vals_sol = load('Solution_Adv_Diff.mat', 'sing_vals').sing_vals;

error = norm(w_sol - w);
error = max(error, norm(sing_vals_sol - sing_vals));
if error ~= 0
    disp('Error in optimal experimental design Adv_Diff example');
end

% save('Solution_Adv_Diff.mat','w','sing_vals')
