%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
clear;
close all;
clc;
addpath(genpath('../../src/'));
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

sen_op = Euclidean_Sensitivity_Operators_Sola(bayes_inv.obj, con);
qn_prec = Quasi_Newton_Preconditioner_Bayesian_Inversion(z_bar, theta_bar, bayes_inv);
sen = Pseudo_Time_Continuation(z_bar, sen_op, qn_prec);

theta_star = 1.0 + 0.2 * (1 - con.x).^2;

%%
rank = 8;
oversampling = 10;
[u, z, lambda, theta] = sen.qn_prec.Compute_Nominal_Hessian(rank, oversampling);
sen.sen_op.current_u = u;
sen.sen_op.current_z = z;
sen.sen_op.current_lambda = lambda;
sen.sen_op.current_theta = theta;

N_fe = 15;
theta_traj = Euclidean_Auxillary_Parameter_Trajectory(N_fe, theta_bar, theta_star);
[z_k_fe, grad_k_fe] = sen.Pseudo_Time_Continuation_Forward_Euler(theta_traj);
num_state_solves_fe = con.num_state_solves;
num_adjoint_solves_fe = con.num_adjoint_solves;

con.num_state_solves = 0;
con.num_adjoint_solves = 0;

%%
rank = 8;
oversampling = 10;
[u, z, lambda, theta] = sen.qn_prec.Compute_Nominal_Hessian(rank, oversampling);
sen.sen_op.current_u = u;
sen.sen_op.current_z = z;
sen.sen_op.current_lambda = lambda;
sen.sen_op.current_theta = theta;

N_me = 10;
theta_traj = Euclidean_Auxillary_Parameter_Trajectory(N_me, theta_bar, theta_star);
[z_k_me, grad_k_me] = sen.Pseudo_Time_Continuation_Modified_Euler(theta_traj);
num_state_solves_me = con.num_state_solves;
num_adjoint_solves_me = con.num_adjoint_solves;

con.num_state_solves = 0;
con.num_adjoint_solves = 0;

%%
grad_norm_nom = norm(sen_op.Euclidean_Gradient(z_bar, theta_bar));
bayes_inv.con.theta_current = theta_star;
[~, z_star] = bayes_inv.Compute_MAP_Point(z_bar);

%%
normalization = sqrt(z_star' * con.M * z_star);
fe_error = sqrt((z_star - z_k_fe(:, end))' * con.M * (z_star - z_k_fe(:, end))) / normalization;
me_error = sqrt((z_star - z_k_me(:, end))' * con.M * (z_star - z_k_me(:, end))) / normalization;

obj_fe = zeros(N_fe + 1, 1);
obj_me = zeros(N_me + 1, 1);
[obj_star, grad_star] = bayes_inv.opt.Jhat(z_star);
for k = 1:(N_fe + 1)
    obj_fe(k) = bayes_inv.opt.Jhat(z_k_fe(:, k));
end
for k = 1:(N_me + 1)
    obj_me(k) = bayes_inv.opt.Jhat(z_k_me(:, k));
end

%%
x = linspace(0, 1, m)';
figure;
hold on;
plot(x, z_star, '-', 'LineWidth', 3);
plot(x, z_k_fe(:, end), '--', 'LineWidth', 3);
plot(x, z_k_me(:, end), ':', 'LineWidth', 3);
legend({'Optimal Solution', 'Forward Euler', 'Modified Euler'}, 'Location', 'northeast', 'FontSize', 20);
xlabel('$x$', 'Interpreter', 'latex');
ylabel('Diffusion Coefficient');
set(gca, 'fontsize', 20);
if write_figure_to_file
    cd(write_path);
    saveas(gca, 'continuation_solution_maps', 'epsc');
    cd(working_path);
end

figure;
hold on;
plot(0:N_fe, obj_star * ones(N_fe + 1, 1), '-', 'LineWidth', 3);
plot(0:N_fe, obj_fe, '--', 'LineWidth', 3);
plot(0:N_me, obj_me, ':', 'LineWidth', 3);
legend({'Optimal Solution', 'Forward Euler', 'Modified Euler'}, 'Location', 'northeast', 'FontSize', 20);
xlabel('Time Step');
ylabel('Objective');
set(gca, 'fontsize', 20);
if write_figure_to_file
    cd(write_path);
    saveas(gca, 'continuation_solution_obj', 'epsc');
    cd(working_path);
end

figure;
semilogy(0:N_fe, grad_norm_nom * ones(N_fe + 1, 1), '-', 'LineWidth', 3);
hold on;
semilogy(0:N_fe, vecnorm(grad_k_fe), '--', 'LineWidth', 3);
semilogy(0:N_me, vecnorm(grad_k_me), '--', 'LineWidth', 3);
legend({'Optimal Solution', 'Forward Euler', 'Modified Euler'}, 'Location', 'best', 'FontSize', 20);
xlabel('Time Step');
ylabel('Gradient Norm');
set(gca, 'fontsize', 20);
if write_figure_to_file
    cd(write_path);
    saveas(gca, 'continuation_solution_grad', 'epsc');
    cd(working_path);
end
