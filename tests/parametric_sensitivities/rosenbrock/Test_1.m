clear;
close all;

execute_fd_check = false;
suppress_plotting = true;

d = 3;
rosenbrock = Rosenbrock(d);

z0 = randn(d, 1);
theta_bar = ones(d - 1, 1);

if execute_fd_check
    rosenbrock.Gradient_FD_Check(z0, theta_bar);
    rosenbrock.Hessian_FD_Check(z0, theta_bar);
    rosenbrock.B_FD_Check(z0, theta_bar);
end

verbose = 'iter';
% z_bar = rosenbrock.Solve_Optimization_Problem(z0, theta_bar, verbose);
z_bar = ones(3, 1);

sen_op = Sensitivity_Operators_Rosenbrock(rosenbrock);
qn_prec = Quasi_Newton_Preconditioner();
sen = Pseudo_Time_Continuation(z_bar, sen_op, qn_prec);

theta_star = 1.2 * ones(d - 1, 1);
evalc('z_star = rosenbrock.Solve_Optimization_Problem(z_bar, theta_star, verbose)');

N_fe = 60;
theta_traj = Euclidean_Auxillary_Parameter_Trajectory(N_fe, theta_bar, theta_star);
evalc('[z_k_fe, grad_k_fe] = sen.Pseudo_Time_Continuation_Forward_Euler(theta_traj)');

N_me = 30;
theta_traj = Euclidean_Auxillary_Parameter_Trajectory(N_me, theta_bar, theta_star);
evalc('[z_k_me, grad_k_me] = sen.Pseudo_Time_Continuation_Modified_Euler(theta_traj)');

error = max(norm(z_star - z_k_fe(:, end)), norm(z_star - z_k_me(:, end)));

if error > 1.e-4
    fprintf(2, '\nparametric_sensitivities/rosenbrock failed.\n');
else
    fprintf(1, '\nparametric_sensitivities/rosenbrock passed.\n');
end

if ~suppress_plotting

    obj_fe = zeros(N_fe + 1, 1);
    obj_me = zeros(N_me + 1, 1);
    [obj_star, grad_star] = rosenbrock.J(z_star, theta_star);
    for k = 1:(N_fe + 1)
        obj_fe(k) = rosenbrock.J(z_k_fe(:, k), theta_star);
    end
    for k = 1:(N_me + 1)
        obj_me(k) = rosenbrock.J(z_k_me(:, k), theta_star);
    end

    figure;
    hold on;
    plot(0:N_fe, obj_fe, '-', 'LineWidth', 3);
    plot(0:N_me, obj_me, '-.', 'LineWidth', 3);
    plot(0:N_fe, obj_star * ones(N_fe + 1, 1), '--', 'LineWidth', 3);
    legend({'Forward Euler', 'Modified Euler', 'Optimal Solution'}, 'Location', 'southeast', 'FontSize', 14);
    xlabel('Iteration', 'FontSize', 18);
    ylabel('Objective', 'FontSize', 18);

    figure;
    semilogy(0:N_fe, vecnorm(grad_k_fe), '-', 'LineWidth', 3);
    hold on;
    semilogy(0:N_me, vecnorm(grad_k_me), '-.', 'LineWidth', 3);
    semilogy(0:N_fe, norm(grad_star) * ones(N_fe + 1, 1), '--', 'LineWidth', 3);
    legend({'Forward Euler', 'Modified Euler', 'Optimal Solution'}, 'Location', 'southeast', 'FontSize', 14);
    xlabel('Iteration', 'FontSize', 18);
    ylabel('Gradient Norm', 'FontSize', 18);

end
