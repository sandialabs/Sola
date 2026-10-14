%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef MD_Continuation_Update < handle

    properties
        post_sampling
        hessian_analysis
        opt_prob_interface
        u_opt
        z_opt
        num_continuation_steps
        r
        discard_cache
        materialize_breve_samples
        use_parallel_sampling
        num_mean_to_sample_steps
        use_mean_time_step_heuristic
        mean_time_step_nonlinearity_tol
        mean_time_step_max_steps
        mean_time_step_options
        last_mean_time_step_prior_diagnostics
        last_mean_time_step_diagnostics
        last_mean_time_grid
        use_sample_time_step_heuristic
        sample_time_step_nonlinearity_tol
        sample_time_step_max_steps
        sample_time_step_options
        last_sample_time_step_prior_diagnostics
        last_sample_time_step_diagnostics
        last_sample_time_grids
    end

    methods

        function this = MD_Continuation_Update(post_sampling, hessian_analysis, num_continuation_steps, discard_cache, materialize_breve_samples, use_parallel_sampling, num_mean_to_sample_steps, use_mean_time_step_heuristic, mean_time_step_nonlinearity_tol, mean_time_step_max_steps, mean_time_step_options, use_sample_time_step_heuristic, sample_time_step_nonlinearity_tol, sample_time_step_max_steps, sample_time_step_options)
            arguments
                post_sampling MD_Posterior_Sampling
                hessian_analysis MD_Hessian_Analysis
                num_continuation_steps (1, 1) {mustBeNumeric}
                discard_cache = true
                materialize_breve_samples (1, 1) logical = false
                use_parallel_sampling (1, 1) logical = false
                num_mean_to_sample_steps (1, 1) {mustBeNumeric} = max(1, ceil(num_continuation_steps / 2))
                use_mean_time_step_heuristic (1, 1) logical = false
                mean_time_step_nonlinearity_tol (1, 1) {mustBeNumeric, mustBeNonnegative} = 0.1
                mean_time_step_max_steps (1, 1) {mustBeNumeric} = max(10 * num_continuation_steps, 20)
                mean_time_step_options struct = struct()
                use_sample_time_step_heuristic (1, 1) logical = false
                sample_time_step_nonlinearity_tol (1, 1) {mustBeNumeric, mustBeNonnegative} = 0.1
                sample_time_step_max_steps (1, 1) {mustBeNumeric} = max(10 * max(1, ceil(num_continuation_steps / 2)), 20)
                sample_time_step_options struct = struct()
            end
            this.post_sampling = post_sampling;
            this.hessian_analysis = hessian_analysis;
            this.opt_prob_interface = hessian_analysis.opt_prob_interface;
            this.u_opt = post_sampling.data_interface.u_opt;
            this.z_opt = post_sampling.data_interface.z_opt;
            if isempty(hessian_analysis.evals)
                this.r = length(this.z_opt);
            else
                this.r = length(hessian_analysis.evals);
            end
            this.num_continuation_steps = num_continuation_steps;
            this.discard_cache = discard_cache;
            this.materialize_breve_samples = materialize_breve_samples;
            this.use_parallel_sampling = use_parallel_sampling;
            this.num_mean_to_sample_steps = num_mean_to_sample_steps;
            this.use_mean_time_step_heuristic = use_mean_time_step_heuristic;
            this.mean_time_step_nonlinearity_tol = mean_time_step_nonlinearity_tol;
            this.mean_time_step_max_steps = mean_time_step_max_steps;
            this.mean_time_step_options = mean_time_step_options;
            this.last_mean_time_step_prior_diagnostics = struct();
            this.last_mean_time_step_diagnostics = {};
            this.last_mean_time_grid = [];
            this.use_sample_time_step_heuristic = use_sample_time_step_heuristic;
            this.sample_time_step_nonlinearity_tol = sample_time_step_nonlinearity_tol;
            this.sample_time_step_max_steps = sample_time_step_max_steps;
            if isempty(fieldnames(sample_time_step_options))
                this.sample_time_step_options = mean_time_step_options;
            else
                this.sample_time_step_options = sample_time_step_options;
            end
            this.last_sample_time_step_prior_diagnostics = {};
            this.last_sample_time_step_diagnostics = {};
            this.last_sample_time_grids = {};
        end

        % ------------------------------------------------------------
        % High-level functions
        % ------------------------------------------------------------

        function [u, z, beta] = Posterior_Update_Mean(this)
            [u, z, beta] = this.Posterior_Update_Core(0);
        end

        function [u_ks, z_ks, beta_ks] = Posterior_Update_Samples(this, use_parallel_sampling)
            if nargin < 2 || isempty(use_parallel_sampling)
                use_parallel_sampling = this.use_parallel_sampling;
            end
            if use_parallel_sampling && ~this.Parallel_Sampling_Available()
                warning('MD_Continuation_Update:ParallelUnavailable', ...
                    ['Parallel sampling requested, but the Parallel Computing Toolbox/DataQueue ', ...
                    'interface is not available. Falling back to serial sampling.']);
                use_parallel_sampling = false;
            end

            num_samples = this.post_sampling.post_data.num_samples;
            u_ks = zeros(length(this.u_opt), num_samples);
            z_ks = zeros(length(this.z_opt), num_samples);
            beta_ks = zeros(this.r, num_samples);

            if use_parallel_sampling
                h = waitbar(0, 'Initializing...', 'Name', 'Processing samples in parallel...');
                tic;
                progress_count = 0;
                q = parallel.pool.DataQueue;
                afterEach(q, @(~) update_parallel_waitbar());

                parfor sample_idx = 1:num_samples
                    [u_k, z_k, beta_k] = this.Posterior_Update_Core(sample_idx);
                    u_ks(:, sample_idx) = u_k;
                    z_ks(:, sample_idx) = z_k;
                    beta_ks(:, sample_idx) = beta_k;
                    send(q, sample_idx);
                end
                close(h);
            else
                h = waitbar(0, 'Initializing...', 'Name', 'Processing...');
                tic;
                for sample_idx = 1:num_samples
                    [u_k, z_k, beta_k] = this.Posterior_Update_Core(sample_idx);
                    u_ks(:, sample_idx) = u_k;
                    z_ks(:, sample_idx) = z_k;
                    beta_ks(:, sample_idx) = beta_k;

                    waitbar(sample_idx/num_samples, h, sprintf('Progress: %d%% | ETA: %.1fs', round(100*sample_idx/num_samples), (toc / sample_idx) * (num_samples - sample_idx)));
                end
                close(h);
            end

            function update_parallel_waitbar()
                progress_count = progress_count + 1;
                waitbar(progress_count/num_samples, h, sprintf('Progress: %d%% | ETA: %.1fs', round(100*progress_count/num_samples), (toc / progress_count) * (num_samples - progress_count)));
            end
        end

        function [u_ks, z_ks, beta_ks, u_bar, z_bar, beta_bar] = Posterior_Update_Samples_From_Mean(this, u_bar, z_bar, beta_bar, use_parallel_sampling)
            if nargin < 2
                u_bar = [];
            end
            if nargin < 3
                z_bar = [];
            end
            if nargin < 4
                beta_bar = [];
            end
            if nargin < 5 || isempty(use_parallel_sampling)
                use_parallel_sampling = this.use_parallel_sampling;
            end
            if use_parallel_sampling && ~this.Parallel_Sampling_Available()
                warning('MD_Continuation_Update:ParallelUnavailable', ...
                    ['Parallel sampling requested, but the Parallel Computing Toolbox/DataQueue ', ...
                    'interface is not available. Falling back to serial sampling.']);
                use_parallel_sampling = false;
            end

            % Continue each posterior sample from the posterior mean update.
            if isempty(beta_bar)
                [u_bar, z_bar, beta_bar] = this.Posterior_Update_Mean();
            else
                if isempty(z_bar)
                    z_bar = this.z_opt + this.hessian_analysis.Apply_V(beta_bar);
                end
                if isempty(u_bar)
                    u_bar = this.opt_prob_interface.State_Solve(z_bar);
                end
            end

            num_samples = this.post_sampling.post_data.num_samples;
            u_ks = zeros(length(this.u_opt), num_samples);
            z_ks = zeros(length(this.z_opt), num_samples);
            beta_ks = zeros(this.r, num_samples);
            prior_diagnostics = cell(num_samples, 1);
            step_diagnostics = cell(num_samples, 1);
            time_grids = cell(num_samples, 1);

            if use_parallel_sampling
                h = waitbar(0, 'Initializing...', 'Name', 'Processing mean-to-sample continuation in parallel...');
                tic;
                progress_count = 0;
                q = parallel.pool.DataQueue;
                afterEach(q, @(~) update_parallel_waitbar());

                parfor sample_idx = 1:num_samples
                    [u_k, z_k, beta_k, prior_diag, time_grid, step_diag] = this.Posterior_Update_Sample_From_Mean_Core(beta_bar, sample_idx);
                    u_ks(:, sample_idx) = u_k;
                    z_ks(:, sample_idx) = z_k;
                    beta_ks(:, sample_idx) = beta_k;
                    prior_diagnostics{sample_idx} = prior_diag;
                    time_grids{sample_idx} = time_grid;
                    step_diagnostics{sample_idx} = step_diag;
                    send(q, sample_idx);
                end
                close(h);
            else
                sen_op = MD_Continuation_Sensitivity_Operators(this.post_sampling, this.hessian_analysis, this.discard_cache, this.materialize_breve_samples);

                h = waitbar(0, 'Initializing...', 'Name', 'Processing mean-to-sample continuation...');
                tic;
                for sample_idx = 1:num_samples
                    [u_k, z_k, beta_k, prior_diag, time_grid, step_diag] = this.Posterior_Update_Sample_From_Mean_Core(beta_bar, sample_idx, sen_op);
                    u_ks(:, sample_idx) = u_k;
                    z_ks(:, sample_idx) = z_k;
                    beta_ks(:, sample_idx) = beta_k;
                    prior_diagnostics{sample_idx} = prior_diag;
                    time_grids{sample_idx} = time_grid;
                    step_diagnostics{sample_idx} = step_diag;

                    waitbar(sample_idx/num_samples, h, sprintf('Progress: %d%% | ETA: %.1fs', round(100*sample_idx/num_samples), (toc / sample_idx) * (num_samples - sample_idx)));
                end
                close(h);
            end

            this.last_sample_time_step_prior_diagnostics = prior_diagnostics;
            this.last_sample_time_grids = time_grids;
            this.last_sample_time_step_diagnostics = step_diagnostics;

            function update_parallel_waitbar()
                progress_count = progress_count + 1;
                waitbar(progress_count/num_samples, h, sprintf('Progress: %d%% | ETA: %.1fs', round(100*progress_count/num_samples), (toc / progress_count) * (num_samples - progress_count)));
            end
        end

        % ------------------------------------------------------------
        % Core continuation driver (shared by mean and samples)
        % ------------------------------------------------------------

        function [u_k, z_k, beta_k] = Posterior_Update_Core(this, sample_idx)

            % NOTE: All continuation classes prefixed with `MD_` update the latent variable beta, 
            % although denoted by parent classes as z. This is not to be confused with the solution
            % space optimization variable z in model discrepancy calibration.
            sen_op = MD_Continuation_Sensitivity_Operators(this.post_sampling, this.hessian_analysis, this.discard_cache, this.materialize_breve_samples);
            qn_prec = MD_Quasi_Newton_Preconditioner(this.hessian_analysis);
            beta_nom = zeros(this.r, 1);
            pt_cont = Pseudo_Time_Continuation(beta_nom, sen_op, qn_prec);
            if sample_idx == 0 && this.use_mean_time_step_heuristic
                % Use a trajectory with N = 1 so the backtracking integrator can
                % pass physical pseudo-times directly to theta_traj.Get_Time().
                theta_traj = MD_Discrepancy_Parameter_Trajectory(1, sample_idx);

                qn_prec.Set_N(2);
                B = sen_op.Apply_B(beta_nom, theta_traj, 0);
                f = pt_cont.Apply_Inverse_Hessian(B, beta_nom, theta_traj, 0);
                beta_dot = -f;
                [dt_prior, this.last_mean_time_step_prior_diagnostics] = ...
                    sen_op.Mean_A_Priori_Time_Step_Heuristic(beta_nom, theta_traj, 0, beta_dot, ...
                                                             this.mean_time_step_nonlinearity_tol, ...
                                                             this.mean_time_step_options);

                backtracking_options = this.mean_time_step_options;
                if ~isfield(backtracking_options, 'max_steps') || isempty(backtracking_options.max_steps)
                    backtracking_options.max_steps = this.mean_time_step_max_steps;
                end
                [beta_k_tmp, ~, this.last_mean_time_grid, this.last_mean_time_step_diagnostics] = ...
                    pt_cont.Pseudo_Time_Continuation_Forward_Euler_Backtracking(theta_traj, dt_prior, ...
                                                                                backtracking_options);
            else
                theta_traj = MD_Discrepancy_Parameter_Trajectory(this.num_continuation_steps, sample_idx);
                beta_k_tmp = pt_cont.Pseudo_Time_Continuation_Forward_Euler(theta_traj);
                if sample_idx == 0
                    this.last_mean_time_step_prior_diagnostics = struct();
                    this.last_mean_time_grid = linspace(0, 1, this.num_continuation_steps + 1);
                    this.last_mean_time_step_diagnostics = {};
                end
            end

            beta_k = beta_k_tmp(:, end);
            u_k = pt_cont.sen_op.current_u;
            z_k = pt_cont.sen_op.current_z;

        end

        function [u_k, z_k, beta_k, prior_diagnostics, time_grid, step_diagnostics] = Posterior_Update_Sample_From_Mean_Core(this, beta_bar, sample_idx, sen_op)
            if nargin < 4 || isempty(sen_op)
                sen_op = MD_Continuation_Sensitivity_Operators(this.post_sampling, this.hessian_analysis, this.discard_cache, this.materialize_breve_samples);
            end

            qn_prec = MD_Quasi_Newton_Preconditioner(this.hessian_analysis);
            pt_cont = Pseudo_Time_Continuation(beta_bar, sen_op, qn_prec);
            prior_diagnostics = struct();
            time_grid = [];
            step_diagnostics = {};

            if this.use_sample_time_step_heuristic
                theta_traj = MD_Discrepancy_Mean_To_Sample_Trajectory(1, sample_idx);

                qn_prec.Set_N(2);
                B = sen_op.Apply_B(beta_bar, theta_traj, 0);
                f = pt_cont.Apply_Inverse_Hessian(B, beta_bar, theta_traj, 0);
                beta_dot = -f;
                [dt_prior, prior_diagnostics] = ...
                    sen_op.Sample_A_Priori_Time_Step_Heuristic(beta_bar, theta_traj, 0, beta_dot, sample_idx, ...
                                                               this.sample_time_step_nonlinearity_tol, ...
                                                               this.sample_time_step_options);

                backtracking_options = this.sample_time_step_options;
                if ~isfield(backtracking_options, 'max_steps') || isempty(backtracking_options.max_steps)
                    backtracking_options.max_steps = this.sample_time_step_max_steps;
                end
                [beta_k_tmp, ~, time_grid, step_diagnostics] = ...
                    pt_cont.Pseudo_Time_Continuation_Forward_Euler_Backtracking(theta_traj, dt_prior, ...
                                                                                backtracking_options);
            else
                theta_traj = MD_Discrepancy_Mean_To_Sample_Trajectory(this.num_mean_to_sample_steps, sample_idx);
                beta_k_tmp = pt_cont.Pseudo_Time_Continuation_Forward_Euler(theta_traj);
                time_grid = linspace(0, 1, this.num_mean_to_sample_steps + 1);
            end

            beta_k = beta_k_tmp(:, end);
            u_k = pt_cont.sen_op.current_u;
            z_k = pt_cont.sen_op.current_z;
        end

        function [tf] = Parallel_Sampling_Available(this) %#ok<MANU>
            tf = exist('parallel.pool.DataQueue', 'class') == 8;
            if tf && exist('license', 'builtin') == 5
                tf = license('test', 'Distrib_Computing_Toolbox');
            end
        end

        % ------------------------------------------------------------
        % Objective helpers (optional; kept for completeness)
        % ------------------------------------------------------------

        function [val, grad] = Parameterized_RS_Objective_beta(this, beta, sample_idx)
            beta = beta(:);
            sen_op = MD_Continuation_Sensitivity_Operators(this.post_sampling, this.hessian_analysis, this.discard_cache, this.materialize_breve_samples);
            disc_ops = sen_op.Get_Continuation_Beta_Discrepancy_Ops(sample_idx);
            z = this.z_opt + this.hessian_analysis.Apply_V(beta);
            u = this.opt_prob_interface.State_Solve(z);
            delta = disc_ops.Eval(beta, 1.0);

            [val, grad_u, grad_z] = this.opt_prob_interface.Objective_Function(u + delta, z);
            beta_grad_z = this.hessian_analysis.Apply_V_Transpose(grad_z);
            z_S_adj = this.opt_prob_interface.Apply_Solution_Operator_z_Jacobian_Transpose(grad_u, z);
            beta_grad_S = this.hessian_analysis.Apply_V_Transpose(z_S_adj);
            beta_grad_D = disc_ops.Apply_Beta_Jacobian_Transpose(grad_u, 1.0);
            grad = beta_grad_z + beta_grad_S + beta_grad_D;

        end

        function [val, grad] = Parameterized_RS_Objective(this, z, sample_idx)

            assert(sample_idx == 0, ...
                ['Parameterized_RS_Objective is only valid for the mean discrepancy. ', ...
                'For posterior samples, use Parameterized_RS_Objective_beta.']);

            sen_op = MD_Continuation_Sensitivity_Operators(this.post_sampling, this.hessian_analysis, this.discard_cache, this.materialize_breve_samples);

            u = this.opt_prob_interface.State_Solve(z);
            delta = sen_op.Discrepancy_Evaluation_Mean(z);

            [val, grad_u, grad_z] = this.opt_prob_interface.Objective_Function(u + delta, z);
            z_tmp1 = sen_op.Apply_Discrepancy_z_Jacobian_Transpose_Mean(grad_u);
            z_tmp2 = this.opt_prob_interface.Apply_Solution_Operator_z_Jacobian_Transpose(grad_u, z);
            grad = grad_z + z_tmp1 + z_tmp2;

        end

    end

end
