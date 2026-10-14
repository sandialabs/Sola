%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Pseudo_Time_Continuation < handle

    properties
        z_bar
        sen_op
        qn_prec
        use_qn_prec
        print_output
        print_iter
    end

    methods

        function this = Pseudo_Time_Continuation(z_bar, sen_op, qn_prec, print_output, print_iter)
            arguments
                z_bar (:, 1) {mustBeNumeric}
                sen_op Sensitivity_Operators
                qn_prec Quasi_Newton_Preconditioner
                print_output logical = false
                print_iter logical = false
            end
            this.z_bar = z_bar;
            this.sen_op = sen_op;
            this.qn_prec = qn_prec;
            this.use_qn_prec = true;
            this.print_output = print_output;
            this.print_iter = print_iter;
        end

        function [z_out] = Apply_Inverse_Hessian(this, z_in, z, theta_traj, time_index)

            tol = 1.e-6;
            max_iter = length(z_in);

            z_out = 0 * z_in;
            r = z_in;
            v = this.qn_prec.Apply_Inverse_Hessian_Approximation(r);
            p = v;
            scalar = r' * v;
            rel_tol = tol * norm(z_in);
            iter = 0;

            while (sqrt(scalar) > rel_tol) && (norm(r) > rel_tol) && (iter < max_iter)
                iter = iter + 1;
                if this.print_iter
                    disp(['Iteration = ', num2str(iter)]);
                    disp(['Sqrt(rho) = ', num2str(sqrt(scalar))]);
                    disp(['Norm(r) = ', num2str(norm(r))]);
                end

                w = this.sen_op.Apply_Hessian(p, z, theta_traj, time_index);

                p_dot_w = w' * p;
                this.qn_prec.Add_Block_Quasi_Newton_Step(p, w, p_dot_w);
                alpha = scalar / p_dot_w;
                z_out = z_out + alpha * p;
                r = r - alpha * w;
                v = this.qn_prec.Apply_Inverse_Hessian_Approximation(r);
                scalar_old = scalar;
                scalar = v' * r;
                p = v + (scalar / scalar_old) * p;
            end

            if this.print_output
                disp(['Total iterations = ', num2str(iter)]);
                disp(['Relative Residual = ', num2str(norm(r) / norm(z_in))]);
            end
        end

        function [z_k, grad_k] = Pseudo_Time_Continuation_Forward_Euler(this, theta_traj)

            N = theta_traj.Get_Number_of_Timesteps();
            this.qn_prec.Set_N(N);

            z_k = zeros(length(this.z_bar), N + 1);
            z_k(:, 1) = this.z_bar;

            grad_k = zeros(length(this.z_bar), N + 1);
            dt = 1 / N;

            grad_k(:, 1) = this.sen_op.Gradient(z_k(:, 1), theta_traj, 0);
            for k = 1:N

                B = this.sen_op.Apply_B(z_k(:, k), theta_traj, k - 1);
                f = this.Apply_Inverse_Hessian(B, z_k(:, k), theta_traj, k - 1);

                z_tmp = z_k(:, k) - dt * f;
                grad_tmp = this.sen_op.Gradient(z_tmp, theta_traj, k);

                if this.use_qn_prec
                    this.qn_prec.Add_Block_Quasi_Newton_Data();
                    s_k = z_tmp - z_k(:, k);
                    y_k = grad_tmp - grad_k(:, k) - dt * B;
                    this.qn_prec.Add_Parametric_Quasi_Newton_Data(s_k, y_k);
                end

                f = this.Apply_Inverse_Hessian(grad_tmp, z_tmp, theta_traj, k);
                z_k(:, k + 1) = z_tmp - f;
                grad_k(:, k + 1) = this.sen_op.Gradient(z_k(:, k + 1), theta_traj, k);

                if this.use_qn_prec
                    this.qn_prec.Add_Block_Quasi_Newton_Data();
                    s_k = z_k(:, k + 1) - z_tmp;
                    y_k = grad_k(:, k + 1) - grad_tmp;
                    this.qn_prec.Add_Parametric_Quasi_Newton_Data(s_k, y_k);
                end

            end
        end

        function [z_k, grad_k, time_grid, step_diagnostics] = Pseudo_Time_Continuation_Forward_Euler_Backtracking(this, theta_traj, initial_dt, options)
            arguments
                this
                theta_traj Auxillary_Parameter_Trajectory
                initial_dt (1, 1) {mustBeNumeric, mustBePositive}
                options struct = struct()
            end

            max_steps = this.Get_Option(options, 'max_steps', 100);
            shrink_factor = this.Get_Option(options, 'shrink_factor', 0.5);
            min_dt = this.Get_Option(options, 'min_dt', 1.e-6);
            max_backtracks = this.Get_Option(options, 'max_backtracks', 20);
            gradient_reduction_factor = this.Get_Option(options, 'gradient_reduction_factor', 0.9);
            gradient_abs_tol = this.Get_Option(options, 'gradient_abs_tol', 1.e-8);
            require_objective_decrease = this.Get_Option(options, 'require_objective_decrease', false);
            objective_armijo_factor = this.Get_Option(options, 'objective_armijo_factor', 0.0);

            if shrink_factor <= 0 || shrink_factor >= 1
                error('Pseudo_Time_Continuation:InvalidShrinkFactor', ...
                      'Backtracking shrink_factor must lie in (0, 1).');
            end
            initial_dt = min(initial_dt, 1.0);
            this.qn_prec.Set_N(2 * max_steps);

            z_k = zeros(length(this.z_bar), max_steps + 1);
            z_k(:, 1) = this.z_bar;

            grad_k = zeros(length(this.z_bar), max_steps + 1);
            time_grid = zeros(1, max_steps + 1);
            step_diagnostics = cell(max_steps, 1);

            trajectory_scale = theta_traj.Get_Number_of_Timesteps();
            t = 0.0;
            grad_k(:, 1) = this.sen_op.Gradient(z_k(:, 1), theta_traj, 0);
            k = 1;

            while t < 1.0 - 10 * eps
                if k > max_steps
                    error('Pseudo_Time_Continuation:BacktrackingMaxSteps', ...
                          'Pseudo-time continuation exceeded max_steps = %d.', max_steps);
                end

                trial_dt = min(initial_dt, 1.0 - t);
                accepted = false;
                num_backtracks = 0;
                last_diagnostics = struct();

                while ~accepted
                    if trial_dt < min_dt
                        error('Pseudo_Time_Continuation:BacktrackingMinStep', ...
                              'Pseudo-time continuation step fell below min_dt = %.3e.', min_dt);
                    end
                    if num_backtracks > max_backtracks
                        error('Pseudo_Time_Continuation:BacktrackingFailure', ...
                              'Pseudo-time continuation exceeded max_backtracks = %d.', max_backtracks);
                    end

                    qn_snapshot = this.qn_prec.Snapshot();
                    time_index = trajectory_scale * t;
                    t_next = t + trial_dt;
                    time_index_next = trajectory_scale * t_next;

                    B = this.sen_op.Apply_B(z_k(:, k), theta_traj, time_index);
                    f = this.Apply_Inverse_Hessian(B, z_k(:, k), theta_traj, time_index);

                    z_pred = z_k(:, k) - trial_dt * f;
                    [grad_pred, val_pred] = this.sen_op.Gradient(z_pred, theta_traj, time_index_next);

                    if this.use_qn_prec
                        this.qn_prec.Add_Block_Quasi_Newton_Data();
                        s_k = z_pred - z_k(:, k);
                        y_k = grad_pred - grad_k(:, k) - trial_dt * B;
                        this.qn_prec.Add_Parametric_Quasi_Newton_Data(s_k, y_k);
                    end

                    f_corr = this.Apply_Inverse_Hessian(grad_pred, z_pred, theta_traj, time_index_next);
                    z_corr = z_pred - f_corr;
                    [grad_corr, val_corr] = this.sen_op.Gradient(z_corr, theta_traj, time_index_next);

                    if this.use_qn_prec
                        this.qn_prec.Add_Block_Quasi_Newton_Data();
                        s_k = z_corr - z_pred;
                        y_k = grad_corr - grad_pred;
                        this.qn_prec.Add_Parametric_Quasi_Newton_Data(s_k, y_k);
                    end

                    grad_pred_norm = norm(grad_pred);
                    grad_corr_norm = norm(grad_corr);
                    gradient_accept = grad_corr_norm <= max(gradient_abs_tol, gradient_reduction_factor * grad_pred_norm);
                    objective_accept = ~isempty(val_pred) && ~isempty(val_corr) && ...
                        isscalar(val_pred) && isscalar(val_corr) && ...
                        isfinite(val_pred) && isfinite(val_corr) && ...
                        (val_corr <= val_pred - objective_armijo_factor * max(1.0, abs(val_pred)));

                    if require_objective_decrease
                        accepted = gradient_accept && objective_accept;
                    else
                        accepted = gradient_accept || objective_accept;
                    end

                    last_diagnostics = struct();
                    last_diagnostics.t_start = t;
                    last_diagnostics.t_end = t_next;
                    last_diagnostics.dt_used = trial_dt;
                    last_diagnostics.num_backtracks = num_backtracks;
                    last_diagnostics.grad_pred_norm = grad_pred_norm;
                    last_diagnostics.grad_corr_norm = grad_corr_norm;
                    last_diagnostics.val_pred = val_pred;
                    last_diagnostics.val_corr = val_corr;
                    last_diagnostics.gradient_accept = gradient_accept;
                    last_diagnostics.objective_accept = objective_accept;
                    last_diagnostics.accepted = accepted;

                    if ~accepted
                        this.qn_prec.Restore(qn_snapshot);
                        trial_dt = shrink_factor * trial_dt;
                        num_backtracks = num_backtracks + 1;
                    end
                end

                z_k(:, k + 1) = z_corr;
                grad_k(:, k + 1) = grad_corr;
                t = t_next;
                time_grid(k + 1) = t;
                step_diagnostics{k} = last_diagnostics;
                k = k + 1;
            end

            z_k = z_k(:, 1:k);
            grad_k = grad_k(:, 1:k);
            time_grid = time_grid(1:k);
            step_diagnostics = step_diagnostics(1:(k - 1));
        end

        function [value] = Get_Option(this, options, name, default_value) %#ok<INUSL>
            if isfield(options, name) && ~isempty(options.(name))
                value = options.(name);
            else
                value = default_value;
            end
        end

        function [z_k, grad_k] = Pseudo_Time_Continuation_Modified_Euler(this, theta_traj)

            N = theta_traj.Get_Number_of_Timesteps();
            this.qn_prec.Set_N(N);

            z_k = zeros(length(this.z_bar), N + 1);
            z_k(:, 1) = this.z_bar;

            grad_k = zeros(length(this.z_bar), N + 1);
            dt = 1 / N;

            grad_k(:, 1) = this.sen_op.Gradient(z_k(:, 1), theta_traj, 0);
            for k = 1:N

                B = this.sen_op.Apply_B(z_k(:, k), theta_traj, k - 1);
                f = this.Apply_Inverse_Hessian(B, z_k(:, k), theta_traj, k - 1);
                z_tmp1 = z_k(:, k) - 0.5 * dt * f;

                grad_tmp1 = this.sen_op.Gradient(z_tmp1, theta_traj, k - .5);

                if this.use_qn_prec
                    this.qn_prec.Add_Block_Quasi_Newton_Data();
                    s_k = z_tmp1 - z_k(:, k);
                    y_k = grad_tmp1 - grad_k(:, k) - 0.5 * dt * B;
                    this.qn_prec.Add_Parametric_Quasi_Newton_Data(s_k, y_k);
                end

                B = this.sen_op.Apply_B(z_tmp1, theta_traj, k - .5);
                f = this.Apply_Inverse_Hessian(B, z_tmp1, theta_traj, k - .5);
                z_tmp2 = z_k(:, k) - dt * f;

                grad_tmp2 = this.sen_op.Gradient(z_tmp2, theta_traj, k);

                if this.use_qn_prec
                    this.qn_prec.Add_Block_Quasi_Newton_Data();
                    s_k = z_tmp2 - z_tmp1;
                    y_k = grad_tmp2 - grad_tmp1 - 0.5 * dt * B;
                    this.qn_prec.Add_Parametric_Quasi_Newton_Data(s_k, y_k);
                end

                f = this.Apply_Inverse_Hessian(grad_tmp2, z_tmp2, theta_traj, k);
                z_k(:, k + 1) = z_tmp2 - f;
                grad_k(:, k + 1) = this.sen_op.Gradient(z_k(:, k + 1), theta_traj, k);

                if this.use_qn_prec
                    this.qn_prec.Add_Block_Quasi_Newton_Data();
                    s_k = z_k(:, k + 1) - z_tmp2;
                    y_k = grad_k(:, k + 1) - grad_tmp2;
                    this.qn_prec.Add_Parametric_Quasi_Newton_Data(s_k, y_k);
                end

            end
        end

    end
end
