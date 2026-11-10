%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = Run_Samples(sen_op, hessian_analysis, z_opt, r, method, N_c, beta_bar, epsilon_nl, bt_options, sample_indices, num_max_newton_steps, newton_gradient_tol, extrap_scale, use_parallel)
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
    % extrap_scale (optional, default 1): for the mean-to-sample methods, continue to
    % mean + extrap_scale * (sample - mean) (MD_Discrepancy_Extrapolated_Trajectory).
    % extrap_scale > 1 makes the problem harder; scoring uses the same endpoint.
    % use_parallel (optional, default true): run the samples with parfor (current
    % pool, else a Threads pool).  false runs them serially in the client.
    if nargin < 14 || isempty(use_parallel)
        use_parallel = true;
    end
    if nargin < 13 || isempty(extrap_scale)
        extrap_scale = 1;
    end
    mts_traj = @(N, idx) MD_Discrepancy_Extrapolated_Trajectory(N, idx, extrap_scale);
    if nargin < 11 || isempty(num_max_newton_steps)
        num_max_newton_steps = 0;
    end
    if nargin < 12 || isempty(newton_gradient_tol)
        newton_gradient_tol = 0;
    end
    num_samples = length(sample_indices);
    o_z = nan(length(z_opt), num_samples);
    o_grad_norm = nan(num_samples, 1);
    o_steps = nan(num_samples, 1);
    o_backtracks = zeros(num_samples, 1);
    o_dt_prior = nan(num_samples, 1);
    % A priori heuristic quantities (adaptive_mts only; NaN otherwise).
    o_A = nan(num_samples, 1);
    o_eta = nan(num_samples, 1);
    o_dt_nl = nan(num_samples, 1);
    o_dt_curv = nan(num_samples, 1);
    o_failures = false(num_samples, 1);
    o_newton_steps = nan(num_samples, 1);
    % Gradient norm after the continuation steps (column 1) and after each
    % Newton step taken (columns 2..); the Newton iterates are prefixes, so with
    % newton_gradient_tol = 0 one run with M steps gives the result for every
    % smaller M.  Columns beyond the steps taken (early stop) stay NaN.
    o_ngn = nan(num_samples, num_max_newton_steps + 1);

    % Samples are independent: run them in parallel.  Use the current pool if there
    % is one, otherwise start a Threads pool (needs the Parallel Computing Toolbox;
    % without it parfor runs serially).  Each worker gets its own copy of sen_op,
    % so its state cache is never shared; the (seeded) breve samplers already
    % stored in sen_op are copied along with it.  out.time is wall-clock time.
    num_workers = 0;                     % parfor (ii = ..., 0) runs serially in the client
    if use_parallel && license('test', 'Distrib_Computing_Toolbox')
        pool = gcp('nocreate');
        if isempty(pool)
            pool = parpool('Threads');
        end
        num_workers = pool.NumWorkers;
    end
    tic;
    parfor (ii = 1:num_samples, num_workers)
        idx = sample_indices(ii);
        warn_state = warning('error', 'MATLAB:nearlySingularMatrix');
        warn_cleanup = onCleanup(@() warning(warn_state));
        try
            qn_prec = MD_Quasi_Newton_Preconditioner(hessian_analysis);
            switch method
                case 'nominal'
                    pt = Pseudo_Time_Continuation(zeros(r, 1), sen_op, qn_prec);
                    theta = MD_Discrepancy_Parameter_Trajectory(N_c, idx);
                    [beta_traj, grad_traj] = pt.Pseudo_Time_Continuation_Forward_Euler(theta, num_max_newton_steps, newton_gradient_tol);
                    cont_steps = N_c;
                    o_steps(ii) = N_c;
                    o_dt_prior(ii) = 1 / N_c;
                case 'fixed_mts'
                    pt = Pseudo_Time_Continuation(beta_bar, sen_op, qn_prec);
                    theta = mts_traj(N_c, idx);
                    [beta_traj, grad_traj] = pt.Pseudo_Time_Continuation_Forward_Euler(theta, num_max_newton_steps, newton_gradient_tol);
                    cont_steps = N_c;
                    o_steps(ii) = N_c;
                    o_dt_prior(ii) = 1 / N_c;
                case 'adaptive_mts'
                    pt = Pseudo_Time_Continuation(beta_bar, sen_op, qn_prec);
                    theta = mts_traj(1, idx);
                    opts = bt_options;
                    opts.num_max_newton_steps = num_max_newton_steps;
                    opts.newton_gradient_tol = newton_gradient_tol;
                    opts.time_step_callback = @(beta, th, ti, beta_dot) ...
                        sen_op.Sample_A_Priori_Time_Step_Heuristic(beta, th, ti, beta_dot, idx, epsilon_nl, bt_options);
                    [beta_traj, grad_traj, ~, step_diags, prior_diag] = pt.Pseudo_Time_Continuation_Forward_Euler_Backtracking(theta, [], opts);
                    cont_steps = length(step_diags);
                    o_dt_prior(ii) = prior_diag.dt;
                    o_A(ii) = prior_diag.A;
                    o_eta(ii) = prior_diag.eta;
                    o_dt_nl(ii) = prior_diag.dt_nonlinearity;
                    o_dt_curv(ii) = prior_diag.dt_curvature;
                    o_steps(ii) = length(step_diags);
                    o_backtracks(ii) = sum(cellfun(@(d) d.num_backtracks, step_diags));
                case 'backtrack_mts'
                    pt = Pseudo_Time_Continuation(beta_bar, sen_op, qn_prec);
                    theta = mts_traj(1, idx);
                    opts = bt_options;
                    opts.num_max_newton_steps = num_max_newton_steps;
                    opts.newton_gradient_tol = newton_gradient_tol;
                    [beta_traj, grad_traj, ~, step_diags] = pt.Pseudo_Time_Continuation_Forward_Euler_Backtracking(theta, 1.0, opts);
                    cont_steps = length(step_diags);
                    o_dt_prior(ii) = 1.0;
                    o_steps(ii) = length(step_diags);
                    o_backtracks(ii) = sum(cellfun(@(d) d.num_backtracks, step_diags));
                otherwise
                    error('Unknown method "%s".', method);
            end
            beta_k = beta_traj(:, end);
            g_newton = sqrt(sum(grad_traj(:, cont_steps + 1:end).^2, 1));
            ngn_row = nan(1, num_max_newton_steps + 1);
            ngn_row(1:numel(g_newton)) = g_newton;
            o_ngn(ii, :) = ngn_row;
            o_newton_steps(ii) = numel(g_newton) - 1;
            o_z(:, ii) = z_opt + hessian_analysis.Apply_V(beta_k);

            % All methods are scored on the same endpoint (t = 1) sample objective.
            if strcmp(method, 'nominal')
                theta_end = MD_Discrepancy_Parameter_Trajectory(1, idx);
            else
                theta_end = mts_traj(1, idx);
            end
            o_grad_norm(ii) = norm(sen_op.Gradient(beta_k, theta_end, 1));
        catch ME
            o_failures(ii) = true;
            fprintf('  %s (N_c = %s, lambda = %g) sample %d failed: %s\n', method, num2str(N_c), extrap_scale, idx, ME.message);
        end
    end
    out.time = toc;
    out.z = o_z;  out.grad_norm = o_grad_norm;  out.steps = o_steps;  out.backtracks = o_backtracks;
    out.dt_prior = o_dt_prior;  out.A = o_A;  out.eta = o_eta;  out.dt_nl = o_dt_nl;  out.dt_curv = o_dt_curv;
    out.failures = o_failures;  out.newton_steps = o_newton_steps;  out.newton_grad_norms = o_ngn;
    out.work = out.steps + out.backtracks + out.newton_steps;
end
