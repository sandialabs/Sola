clear;
close all;
clc;
addpath('../../src/');
rng(154);

write_figure_to_file = false;
working_path = pwd;
write_path = '~/Desktop';

m = 100;
theta_bar = load('Optimization_Results.mat', 'theta').theta;
con = Adv_Diff_Constraint(theta_bar);
Generate_Obs_Data(con);
likelihood = Adv_Diff_Likelihood_Model(m);
prior = Adv_Diff_Prior_Model(con);
bayes_inv = Bayesian_Inversion(likelihood, prior, con);
bayes_inv.opt.fun_tol = 1.e-10;
bayes_inv.opt.opt_tol = 1.e-10;
z_bar = load('Optimization_Results.mat', 'z_opt').z_opt;

sen_op = Sensitivity_Operators_Sabl(bayes_inv.obj, con);
sen = Pseudo_Time_Continuation_Bayesian_Inversion(z_bar, theta_bar, sen_op, bayes_inv);

theta_star = 1.0 + 0.2 * (1 - con.x).^2;

[z_k_linear_approx, grad_k_linear_approx] = sen.Pseudo_Time_Continuation_Forward_Euler(theta_star, 1);
num_state_solves_linear_approx = con.num_state_solves;
num_adjoint_solves_linear_approx = con.num_adjoint_solves;
linear_approx_cost = num_state_solves_linear_approx + num_adjoint_solves_linear_approx;

%%
N_range = (6:2:30)';
fe_error = 0 * N_range;
fe_cost = 0 * N_range;
fe_bfgs_error = 0 * N_range;
fe_bfgs_cost = 0 * N_range;
me_error = 0 * N_range;
me_cost = 0 * N_range;
me_bfgs_error = 0 * N_range;
me_bfgs_cost = 0 * N_range;

bayes_inv.con.theta_current = theta_star;
[~, z_star] = bayes_inv.Compute_MAP_Point(z_bar);
x = linspace(0, 1, m)';

normalization = sqrt(z_star' * con.M * z_star);
linear_approx_error = sqrt((z_star - z_k_linear_approx(:, end))' * con.M * (z_star - z_k_linear_approx(:, end))) / normalization;

bayes_inv.con.theta_current = theta_bar;

%% Without BFGS
for k = 1:length(N_range)

    %% FE without BFGS
    con.num_state_solves = 0;
    con.num_adjoint_solves = 0;
    sen.use_bfgs_prec = false;
    N_fe = N_range(k);
    [z_k_fe, grad_k_fe] = sen.Pseudo_Time_Continuation_Forward_Euler(theta_star, N_fe);
    num_state_solves_fe = con.num_state_solves;
    num_adjoint_solves_fe = con.num_adjoint_solves;
    fe_error(k) = sqrt((z_star - z_k_fe(:, end))' * con.M * (z_star - z_k_fe(:, end))) / normalization;
    fe_cost(k) = num_state_solves_fe + num_adjoint_solves_fe;

    %% ME without BFGS
    con.num_state_solves = 0;
    con.num_adjoint_solves = 0;
    sen.use_bfgs_prec = false;
    N_me =  N_range(k) / 2;
    [z_k_me, grad_k_me] = sen.Pseudo_Time_Continuation_Modified_Euler(theta_star, N_me);
    num_state_solves_me = con.num_state_solves;
    num_adjoint_solves_me = con.num_adjoint_solves;
    me_error(k) = sqrt((z_star - z_k_me(:, end))' * con.M * (z_star - z_k_me(:, end))) / normalization;
    me_cost(k) = num_state_solves_me + num_adjoint_solves_me;

end

%% With BFGS
for k = 1:length(N_range)
    %% FE with BFGS
    con.num_state_solves = 0;
    con.num_adjoint_solves = 0;
    sen.use_bfgs_prec = true;
    rank = 8;
    oversampling = 10;
    sen.Compute_Nominal_Hessian(rank, oversampling);
    N_fe =  N_range(k);
    [z_k_fe, grad_k_fe] = sen.Pseudo_Time_Continuation_Forward_Euler(theta_star, N_fe);
    num_state_solves_fe = con.num_state_solves;
    num_adjoint_solves_fe = con.num_adjoint_solves;
    fe_bfgs_error(k) = sqrt((z_star - z_k_fe(:, end))' * con.M * (z_star - z_k_fe(:, end))) / normalization;
    fe_bfgs_cost(k) = num_state_solves_fe + num_adjoint_solves_fe;

    %% ME with BFGS
    con.num_state_solves = 0;
    con.num_adjoint_solves = 0;
    sen.use_bfgs_prec = true;
    rank = 8;
    oversampling = 10;
    sen.Compute_Nominal_Hessian(rank, oversampling);
    N_me = N_range(k) / 2;
    [z_k_me, grad_k_me] = sen.Pseudo_Time_Continuation_Modified_Euler(theta_star, N_me);
    num_state_solves_me = con.num_state_solves;
    num_adjoint_solves_me = con.num_adjoint_solves;
    me_bfgs_error(k) = sqrt((z_star - z_k_me(:, end))' * con.M * (z_star - z_k_me(:, end))) / normalization;
    me_bfgs_cost(k) = num_state_solves_me + num_adjoint_solves_me;

end

%%
figure;
hold on;
plot(fe_cost, fe_error, '--', 'LineWidth', 3, 'color', "#0072BD");
plot(fe_bfgs_cost, fe_bfgs_error, '-', 'LineWidth', 3, 'color', "#0072BD");
plot(me_cost, me_error, '--', 'LineWidth', 3, 'color', "#D95319");
plot(me_bfgs_cost, me_bfgs_error, '-', 'LineWidth', 3, 'color', "#D95319");
legend({'FE', 'FE+BFGS', 'ME', 'ME+BFGS'});
xlim([170, 870]);
xlabel('Cost (# of PDE solves)');
ylabel('Relative Error');
set(gca, 'fontsize', 20);
if write_figure_to_file
    cd(write_path);
    saveas(gca, 'cost_error_comparison', 'epsc');
    cd(working_path);
end
