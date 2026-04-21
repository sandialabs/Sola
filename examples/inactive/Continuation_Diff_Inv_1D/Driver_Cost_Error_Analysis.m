%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
close all;
clc;
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

theta_star = 1.0 + 0.2 * (1 - con.x).^2;

sen_op = Euclidean_Sensitivity_Operators_Sabl(bayes_inv.obj, con);
qn_prec = Quasi_Newton_Preconditioner_Bayesian_Inversion(z_bar, theta_bar, bayes_inv);
sen = Pseudo_Time_Continuation(z_bar, sen_op, qn_prec);
sen.use_qn_prec = false;
theta_traj = Euclidean_Auxillary_Parameter_Trajectory(1, theta_bar, theta_star);
[z_k_linear_approx, grad_k_linear_approx] = sen.Pseudo_Time_Continuation_Forward_Euler(theta_traj);

num_state_solves_linear_approx = con.num_state_solves;
num_adjoint_solves_linear_approx = con.num_adjoint_solves;
linear_approx_cost = num_state_solves_linear_approx + num_adjoint_solves_linear_approx;

%%
N_range = (6:2:16)';
fe_error = 0 * N_range;
fe_cost = 0 * N_range;
fe_param_qn_error = 0 * N_range;
fe_param_qn_cost = 0 * N_range;
fe_qn_error = 0 * N_range;
fe_qn_cost = 0 * N_range;
me_error = 0 * N_range;
me_cost = 0 * N_range;
me_param_qn_error = 0 * N_range;
me_param_qn_cost = 0 * N_range;
me_qn_error = 0 * N_range;
me_qn_cost = 0 * N_range;

bayes_inv.con.theta_current = theta_star;
[~, z_star] = bayes_inv.Compute_MAP_Point(z_bar);
x = linspace(0, 1, m)';

normalization = sqrt(z_star' * con.M * z_star);
linear_approx_error = sqrt((z_star - z_k_linear_approx(:, end))' * con.M * (z_star - z_k_linear_approx(:, end))) / normalization;

bayes_inv.con.theta_current = theta_bar;

%%
for k = 1:length(N_range)

    %% FE without QN
    disp(['Working on FE without QN with k = ', num2str(k)]);
    con.num_state_solves = 0;
    con.num_adjoint_solves = 0;

    sen_op = Euclidean_Sensitivity_Operators_Sabl(bayes_inv.obj, con);
    qn_prec = Quasi_Newton_Preconditioner_Bayesian_Inversion(z_bar, theta_bar, bayes_inv);
    sen = Pseudo_Time_Continuation(z_bar, sen_op, qn_prec);
    sen.use_qn_prec = false;

    N_fe = N_range(k);
    theta_traj = Euclidean_Auxillary_Parameter_Trajectory(N_fe, theta_bar, theta_star);
    [z_k_fe, grad_k_fe] = sen.Pseudo_Time_Continuation_Forward_Euler(theta_traj);

    num_state_solves_fe = con.num_state_solves;
    num_adjoint_solves_fe = con.num_adjoint_solves;
    fe_error(k) = sqrt((z_star - z_k_fe(:, end))' * con.M * (z_star - z_k_fe(:, end))) / normalization;
    fe_cost(k) = num_state_solves_fe + num_adjoint_solves_fe;

    %% ME without QN
    disp(['Working on ME without QN with k = ', num2str(k)]);
    con.num_state_solves = 0;
    con.num_adjoint_solves = 0;

    sen_op = Euclidean_Sensitivity_Operators_Sabl(bayes_inv.obj, con);
    qn_prec = Quasi_Newton_Preconditioner_Bayesian_Inversion(z_bar, theta_bar, bayes_inv);
    sen = Pseudo_Time_Continuation(z_bar, sen_op, qn_prec);
    sen.use_qn_prec = false;

    N_me =  N_range(k) / 2;
    theta_traj = Euclidean_Auxillary_Parameter_Trajectory(N_me, theta_bar, theta_star);
    [z_k_me, grad_k_me] = sen.Pseudo_Time_Continuation_Modified_Euler(theta_traj);

    num_state_solves_me = con.num_state_solves;
    num_adjoint_solves_me = con.num_adjoint_solves;
    me_error(k) = sqrt((z_star - z_k_me(:, end))' * con.M * (z_star - z_k_me(:, end))) / normalization;
    me_cost(k) = num_state_solves_me + num_adjoint_solves_me;

    %% FE with only parametric QN
    disp(['Working on FE with PQN with k = ', num2str(k)]);
    con.num_state_solves = 0;
    con.num_adjoint_solves = 0;

    sen_op = Euclidean_Sensitivity_Operators_Sabl(bayes_inv.obj, con);
    qn_prec = Quasi_Newton_Preconditioner_Bayesian_Inversion(z_bar, theta_bar, bayes_inv);
    qn_prec.max_size = 0;
    sen = Pseudo_Time_Continuation(z_bar, sen_op, qn_prec);

    rank = 8;
    oversampling = 10;
    [u, z, lambda, theta] = sen.qn_prec.Compute_Nominal_Hessian(rank, oversampling);
    sen.sen_op.current_u = u;
    sen.sen_op.current_z = z;
    sen.sen_op.current_lambda = lambda;
    sen.sen_op.current_theta = theta;

    N_fe =  N_range(k);
    theta_traj = Euclidean_Auxillary_Parameter_Trajectory(N_fe, theta_bar, theta_star);
    [z_k_fe, grad_k_fe] = sen.Pseudo_Time_Continuation_Forward_Euler(theta_traj);

    num_state_solves_fe = con.num_state_solves;
    num_adjoint_solves_fe = con.num_adjoint_solves;
    fe_param_qn_error(k) = sqrt((z_star - z_k_fe(:, end))' * con.M * (z_star - z_k_fe(:, end))) / normalization;
    fe_param_qn_cost(k) = num_state_solves_fe + num_adjoint_solves_fe;

    %% ME with only parametric QN
    disp(['Working on ME with PQN with k = ', num2str(k)]);
    con.num_state_solves = 0;
    con.num_adjoint_solves = 0;

    sen_op = Euclidean_Sensitivity_Operators_Sabl(bayes_inv.obj, con);
    qn_prec = Quasi_Newton_Preconditioner_Bayesian_Inversion(z_bar, theta_bar, bayes_inv);
    qn_prec.max_size = 0;
    sen = Pseudo_Time_Continuation(z_bar, sen_op, qn_prec);

    rank = 8;
    oversampling = 10;
    [u, z, lambda, theta] = sen.qn_prec.Compute_Nominal_Hessian(rank, oversampling);
    sen.sen_op.current_u = u;
    sen.sen_op.current_z = z;
    sen.sen_op.current_lambda = lambda;
    sen.sen_op.current_theta = theta;

    N_me =  N_range(k) / 2;
    theta_traj = Euclidean_Auxillary_Parameter_Trajectory(N_me, theta_bar, theta_star);
    [z_k_me, grad_k_me] = sen.Pseudo_Time_Continuation_Modified_Euler(theta_traj);

    num_state_solves_me = con.num_state_solves;
    num_adjoint_solves_me = con.num_adjoint_solves;
    me_param_qn_error(k) = sqrt((z_star - z_k_me(:, end))' * con.M * (z_star - z_k_me(:, end))) / normalization;
    me_param_qn_cost(k) = num_state_solves_me + num_adjoint_solves_me;

    %% FE with QN
    disp(['Working on FE with QN with k = ', num2str(k)]);
    con.num_state_solves = 0;
    con.num_adjoint_solves = 0;

    sen_op = Euclidean_Sensitivity_Operators_Sabl(bayes_inv.obj, con);
    qn_prec = Quasi_Newton_Preconditioner_Bayesian_Inversion(z_bar, theta_bar, bayes_inv);
    sen = Pseudo_Time_Continuation(z_bar, sen_op, qn_prec);

    rank = 8;
    oversampling = 10;
    [u, z, lambda, theta] = sen.qn_prec.Compute_Nominal_Hessian(rank, oversampling);
    sen.sen_op.current_u = u;
    sen.sen_op.current_z = z;
    sen.sen_op.current_lambda = lambda;
    sen.sen_op.current_theta = theta;

    N_fe =  N_range(k);
    theta_traj = Euclidean_Auxillary_Parameter_Trajectory(N_fe, theta_bar, theta_star);
    [z_k_fe, grad_k_fe] = sen.Pseudo_Time_Continuation_Forward_Euler(theta_traj);

    num_state_solves_fe = con.num_state_solves;
    num_adjoint_solves_fe = con.num_adjoint_solves;
    fe_qn_error(k) = sqrt((z_star - z_k_fe(:, end))' * con.M * (z_star - z_k_fe(:, end))) / normalization;
    fe_qn_cost(k) = num_state_solves_fe + num_adjoint_solves_fe;

    %% ME with QN
    disp(['Working on ME with QN with k = ', num2str(k)]);
    con.num_state_solves = 0;
    con.num_adjoint_solves = 0;

    sen_op = Euclidean_Sensitivity_Operators_Sabl(bayes_inv.obj, con);
    qn_prec = Quasi_Newton_Preconditioner_Bayesian_Inversion(z_bar, theta_bar, bayes_inv);
    sen = Pseudo_Time_Continuation(z_bar, sen_op, qn_prec);

    rank = 8;
    oversampling = 10;
    [u, z, lambda, theta] = sen.qn_prec.Compute_Nominal_Hessian(rank, oversampling);
    sen.sen_op.current_u = u;
    sen.sen_op.current_z = z;
    sen.sen_op.current_lambda = lambda;
    sen.sen_op.current_theta = theta;

    N_me = N_range(k) / 2;
    theta_traj = Euclidean_Auxillary_Parameter_Trajectory(N_me, theta_bar, theta_star);
    [z_k_me, grad_k_me] = sen.Pseudo_Time_Continuation_Modified_Euler(theta_traj);

    num_state_solves_me = con.num_state_solves;
    num_adjoint_solves_me = con.num_adjoint_solves;
    me_qn_error(k) = sqrt((z_star - z_k_me(:, end))' * con.M * (z_star - z_k_me(:, end))) / normalization;
    me_qn_cost(k) = num_state_solves_me + num_adjoint_solves_me;

end

%%
figure;
hold on;
plot(fe_cost, fe_error, '--', 'LineWidth', 3, 'color', "#0072BD");
plot(fe_param_qn_cost, fe_param_qn_error, ':', 'LineWidth', 3, 'color', "#0072BD");
plot(fe_qn_cost, fe_qn_error, '-', 'LineWidth', 3, 'color', "#0072BD");
plot(me_cost, me_error, '--', 'LineWidth', 3, 'color', "#D95319");
plot(me_param_qn_cost, me_param_qn_error, ':', 'LineWidth', 3, 'color', "#D95319");
plot(me_qn_cost, me_qn_error, '-', 'LineWidth', 3, 'color', "#D95319");
legend({'FE', 'FE+PQN', 'FE+FQN', 'ME', 'ME+PQN', 'ME+FQN'});
xlim([170, 870]);
xlabel('Cost (# of PDE solves)');
ylabel('Relative Error');
set(gca, 'fontsize', 20);
if write_figure_to_file
    cd(write_path);
    saveas(gca, 'cost_error_comparison', 'epsc');
    cd(working_path);
end

save('Cost_Error_Analysis.mat');
