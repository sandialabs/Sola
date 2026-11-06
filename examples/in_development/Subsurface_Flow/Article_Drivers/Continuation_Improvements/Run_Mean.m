%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = Run_Mean(sen_op, hessian_analysis, z_opt, r, num_steps, epsilon_nl, bt_options, num_newton_steps)
    % Posterior-mean continuation: fixed uniform steps if epsilon_nl is
    % empty, otherwise the a priori heuristic with backtracking.
    % num_newton_steps (optional, default 0): extra Newton steps at t = 1 after
    % the continuation steps (both modes); they count toward out.work.
    if nargin < 8 || isempty(num_newton_steps)
        num_newton_steps = 0;
    end
    qn_prec = MD_Quasi_Newton_Preconditioner(hessian_analysis);
    pt = Pseudo_Time_Continuation(zeros(r, 1), sen_op, qn_prec);
    tic;
    if isempty(epsilon_nl)
        theta = MD_Discrepancy_Parameter_Trajectory(num_steps, 0);
        beta_traj = pt.Pseudo_Time_Continuation_Forward_Euler(theta, num_newton_steps);
        out.steps = num_steps;
        out.backtracks = 0;
        out.time_grid = linspace(0, 1, num_steps + 1);
        out.dt_prior = 1 / num_steps;
    else
        theta = MD_Discrepancy_Parameter_Trajectory(1, 0);
        bt_options.num_newton_steps = num_newton_steps;
        bt_options.time_step_callback = @(beta, th, ti, beta_dot) ...
            sen_op.Mean_A_Priori_Time_Step_Heuristic(beta, th, ti, beta_dot, epsilon_nl, bt_options);
        [beta_traj, ~, out.time_grid, step_diags, prior_diag] = pt.Pseudo_Time_Continuation_Forward_Euler_Backtracking(theta, [], bt_options);
        out.dt_prior = prior_diag.dt;
        out.steps = length(step_diags);
        out.backtracks = sum(cellfun(@(d) d.num_backtracks, step_diags));
    end
    out.time = toc;
    out.work = out.steps + out.backtracks + num_newton_steps;
    out.beta = beta_traj(:, end);
    out.z = z_opt + hessian_analysis.Apply_V(out.beta);
    grad = sen_op.Gradient(out.beta, theta, theta.Get_Number_of_Timesteps());
    out.grad_norm = norm(grad);
end
