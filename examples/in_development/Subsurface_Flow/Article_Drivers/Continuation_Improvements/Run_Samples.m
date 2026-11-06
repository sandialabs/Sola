%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = Run_Samples(sen_op, hessian_analysis, z_opt, r, method, N_c, beta_bar, epsilon_nl, bt_options, sample_indices, num_max_newton_steps, newton_gradient_tol)
    % method: 'nominal'      fixed steps from beta = 0
    %         'fixed_mts'    fixed steps from beta_bar
    %         'adaptive_mts' heuristic + backtracking from beta_bar
    %         'backtrack_mts' backtracking only: start with a full step (dt = 1) from
    %                         beta_bar and shrink only when a step is rejected
    % num_max_newton_steps (optional, default 0): maximum number of Newton steps
    % at t = 1 after the continuation steps (all methods).
    % newton_gradient_tol (optional, default 0): stop the Newton steps (per
    % sample) once ||g_beta|| <= tol (0: always take num_max_newton_steps).
    % The Newton steps actually taken are out.newton_steps and count toward out.work.
    if nargin < 11 || isempty(num_max_newton_steps)
        num_max_newton_steps = 0;
    end
    if nargin < 12 || isempty(newton_gradient_tol)
        newton_gradient_tol = 0;
    end
    num_samples = length(sample_indices);
    out.z = nan(length(z_opt), num_samples);
    out.grad_norm = nan(num_samples, 1);
    out.steps = nan(num_samples, 1);
    out.backtracks = zeros(num_samples, 1);
    out.dt_prior = nan(num_samples, 1);
    % A priori heuristic quantities (adaptive_mts only; NaN otherwise).
    out.A = nan(num_samples, 1);
    out.eta = nan(num_samples, 1);
    out.dt_nl = nan(num_samples, 1);
    out.dt_curv = nan(num_samples, 1);
    out.failures = false(num_samples, 1);
    out.newton_steps = nan(num_samples, 1);
    % Gradient norm after the continuation steps (column 1) and after each
    % Newton step taken (columns 2..); the Newton iterates are prefixes, so with
    % newton_gradient_tol = 0 one run with M steps gives the result for every
    % smaller M.  Columns beyond the steps taken (early stop) stay NaN.
    out.newton_grad_norms = nan(num_samples, num_max_newton_steps + 1);

    tic;
    for ii = 1:num_samples
        idx = sample_indices(ii);
        try
            qn_prec = MD_Quasi_Newton_Preconditioner(hessian_analysis);
            switch method
                case 'nominal'
                    pt = Pseudo_Time_Continuation(zeros(r, 1), sen_op, qn_prec);
                    theta = MD_Discrepancy_Parameter_Trajectory(N_c, idx);
                    [beta_traj, grad_traj] = pt.Pseudo_Time_Continuation_Forward_Euler(theta, num_max_newton_steps, newton_gradient_tol);
                    cont_steps = N_c;
                    out.steps(ii) = N_c;
                    out.dt_prior(ii) = 1 / N_c;
                case 'fixed_mts'
                    pt = Pseudo_Time_Continuation(beta_bar, sen_op, qn_prec);
                    theta = MD_Discrepancy_Mean_To_Sample_Trajectory(N_c, idx);
                    [beta_traj, grad_traj] = pt.Pseudo_Time_Continuation_Forward_Euler(theta, num_max_newton_steps, newton_gradient_tol);
                    cont_steps = N_c;
                    out.steps(ii) = N_c;
                    out.dt_prior(ii) = 1 / N_c;
                case 'adaptive_mts'
                    pt = Pseudo_Time_Continuation(beta_bar, sen_op, qn_prec);
                    theta = MD_Discrepancy_Mean_To_Sample_Trajectory(1, idx);
                    opts = bt_options;
                    opts.num_max_newton_steps = num_max_newton_steps;
                    opts.newton_gradient_tol = newton_gradient_tol;
                    opts.time_step_callback = @(beta, th, ti, beta_dot) ...
                        sen_op.Sample_A_Priori_Time_Step_Heuristic(beta, th, ti, beta_dot, idx, epsilon_nl, bt_options);
                    [beta_traj, grad_traj, ~, step_diags, prior_diag] = pt.Pseudo_Time_Continuation_Forward_Euler_Backtracking(theta, [], opts);
                    cont_steps = length(step_diags);
                    out.dt_prior(ii) = prior_diag.dt;
                    out.A(ii) = prior_diag.A;
                    out.eta(ii) = prior_diag.eta;
                    out.dt_nl(ii) = prior_diag.dt_nonlinearity;
                    out.dt_curv(ii) = prior_diag.dt_curvature;
                    out.steps(ii) = length(step_diags);
                    out.backtracks(ii) = sum(cellfun(@(d) d.num_backtracks, step_diags));
                case 'backtrack_mts'
                    pt = Pseudo_Time_Continuation(beta_bar, sen_op, qn_prec);
                    theta = MD_Discrepancy_Mean_To_Sample_Trajectory(1, idx);
                    opts = bt_options;
                    opts.num_max_newton_steps = num_max_newton_steps;
                    opts.newton_gradient_tol = newton_gradient_tol;
                    [beta_traj, grad_traj, ~, step_diags] = pt.Pseudo_Time_Continuation_Forward_Euler_Backtracking(theta, 1.0, opts);
                    cont_steps = length(step_diags);
                    out.dt_prior(ii) = 1.0;
                    out.steps(ii) = length(step_diags);
                    out.backtracks(ii) = sum(cellfun(@(d) d.num_backtracks, step_diags));
                otherwise
                    error('Unknown method "%s".', method);
            end
            beta_k = beta_traj(:, end);
            g_newton = sqrt(sum(grad_traj(:, cont_steps + 1:end).^2, 1));
            out.newton_grad_norms(ii, 1:numel(g_newton)) = g_newton;
            out.newton_steps(ii) = numel(g_newton) - 1;
            out.z(:, ii) = z_opt + hessian_analysis.Apply_V(beta_k);

            % All methods are scored on the same endpoint (t = 1) sample objective.
            theta_end = MD_Discrepancy_Parameter_Trajectory(1, idx);
            out.grad_norm(ii) = norm(sen_op.Gradient(beta_k, theta_end, 1));
        catch ME
            out.failures(ii) = true;
            fprintf('  %s sample %d failed: %s\n', method, idx, ME.message);
        end
    end
    out.time = toc;
    out.work = out.steps + out.backtracks + out.newton_steps;
end
