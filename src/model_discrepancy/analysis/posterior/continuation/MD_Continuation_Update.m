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
        r
        options
        last_mean_time_step_prior_diagnostics
        last_mean_time_step_diagnostics
        last_mean_time_grid
        last_sample_time_step_prior_diagnostics
        last_sample_time_step_diagnostics
        last_sample_time_grids
    end

    methods

        function this = MD_Continuation_Update(post_sampling, hessian_analysis)
            arguments
                post_sampling MD_Posterior_Sampling
                hessian_analysis MD_Hessian_Analysis
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

            this.options = this.Default_Options();

            this.last_mean_time_step_prior_diagnostics = struct();
            this.last_mean_time_step_diagnostics = {};
            this.last_mean_time_grid = [];
            this.last_sample_time_step_prior_diagnostics = {};
            this.last_sample_time_step_diagnostics = {};
            this.last_sample_time_grids = {};
        end

        function [options] = Default_Options(this) %#ok<MANU>
            options = struct();

            options.mean = struct();
            options.mean.num_steps = 3;
            options.mean.time_step = struct();
            options.mean.time_step.use = false;
            options.mean.time_step.epsilon_nl = 0.1;
            options.mean.time_step.options = struct();
            options.mean.time_step.options.max_steps = 30;
            options.mean.time_step.options.shrink_factor = 0.5;
            options.mean.time_step.options.min_dt = 1.e-6;
            options.mean.time_step.options.max_backtracks = 5;
            options.mean.time_step.options.gradient_reduction_factor = 0.9;
            options.mean.time_step.options.gradient_abs_tol = 1.e-8;
            options.mean.time_step.options.require_objective_decrease = false;
            options.mean.time_step.options.objective_armijo_factor = 0.0;
            options.mean.time_step.options.fd_rel_step = 1.e-4;
            options.mean.time_step.options.min_denominator = 1.e-14;
            options.mean.time_step.options.max_dt = 1.0;
            options.mean.time_step.options.epsilon_curv = 1.0;

            options.samples = struct();
            options.samples.num_steps = 3;
            options.samples.num_mean_to_sample_steps = 2;
            options.samples.discard_cache = true;
            options.samples.materialize_breve_samples = false;
            options.samples.use_parallel_sampling = false;
            options.samples.time_step = struct();
            options.samples.time_step.use = false;
            options.samples.time_step.epsilon_nl = 0.1;
            options.samples.time_step.options = struct();
            options.samples.time_step.options.max_steps = 15;
            options.samples.time_step.options.shrink_factor = 0.5;
            options.samples.time_step.options.min_dt = 1.e-6;
            options.samples.time_step.options.max_backtracks = 5;
            options.samples.time_step.options.gradient_reduction_factor = 0.9;
            options.samples.time_step.options.gradient_abs_tol = 1.e-8;
            options.samples.time_step.options.require_objective_decrease = false;
            options.samples.time_step.options.objective_armijo_factor = 0.0;
            options.samples.time_step.options.fd_rel_step = 1.e-4;
            options.samples.time_step.options.min_denominator = 1.e-14;
            options.samples.time_step.options.max_dt = 1.0;
            options.samples.time_step.options.epsilon_curv = 1.0;
        end

        function [options] = Merge_Options(this, defaults, user_options) %#ok<INUSL>
            options = defaults;
            if isempty(user_options)
                return;
            end

            user_fields = fieldnames(user_options);
            for i = 1:length(user_fields)
                field_name = user_fields{i};
                if isfield(options, field_name) && isstruct(options.(field_name)) && isstruct(user_options.(field_name))
                    options.(field_name) = this.Merge_Options(options.(field_name), user_options.(field_name));
                else
                    options.(field_name) = user_options.(field_name);
                end
            end
        end

        % ------------------------------------------------------------
        % High-level functions
        % ------------------------------------------------------------

        function [u, z, beta] = Posterior_Update_Mean(this, mean_options)
            if nargin >= 2 && ~isempty(mean_options)
                this.options.mean = this.Merge_Options(this.options.mean, mean_options);
            end
            [u, z, beta] = this.Posterior_Update_Core(0);
        end

        function [u_ks, z_ks, beta_ks] = Posterior_Update_Samples(this, sample_options)
            if nargin >= 2 && ~isempty(sample_options)
                this.options.samples = this.Merge_Options(this.options.samples, sample_options);
            end
            use_parallel_sampling = this.options.samples.use_parallel_sampling;
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

        function [u_ks, z_ks, beta_ks, u_bar, z_bar, beta_bar] = Posterior_Update_Samples_From_Mean(this, u_bar, z_bar, beta_bar, sample_options)
            if nargin < 2
                u_bar = [];
            end
            if nargin < 3
                z_bar = [];
            end
            if nargin < 4
                beta_bar = [];
            end
            if nargin >= 5 && ~isempty(sample_options)
                this.options.samples = this.Merge_Options(this.options.samples, sample_options);
            end

            use_parallel_sampling = this.options.samples.use_parallel_sampling;
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
                sen_op = MD_Continuation_Sensitivity_Operators(this.post_sampling, this.hessian_analysis, this.options.samples.discard_cache, this.options.samples.materialize_breve_samples);

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
            if sample_idx == 0
                continuation_options = this.options.mean;
                sen_op = MD_Continuation_Sensitivity_Operators(this.post_sampling, this.hessian_analysis);
            else
                continuation_options = this.options.samples;
                sen_op = MD_Continuation_Sensitivity_Operators(this.post_sampling, this.hessian_analysis, ...
                                                               continuation_options.discard_cache, ...
                                                               continuation_options.materialize_breve_samples);
            end
            qn_prec = MD_Quasi_Newton_Preconditioner(this.hessian_analysis);
            beta_nom = zeros(this.r, 1);
            pt_cont = Pseudo_Time_Continuation(beta_nom, sen_op, qn_prec);

            if sample_idx == 0 && continuation_options.time_step.use
                % Use a trajectory with N = 1 so the backtracking integrator can
                % pass physical pseudo-times directly to theta_traj.Get_Time().
                theta_traj = MD_Discrepancy_Parameter_Trajectory(1, sample_idx);

                backtracking_options = continuation_options.time_step.options;
                backtracking_options.time_step_callback = ...
                    @(beta, theta, time_index, beta_dot) sen_op.Mean_A_Priori_Time_Step_Heuristic(beta, theta, time_index, beta_dot, ...
                                                                                                  continuation_options.time_step.epsilon_nl, ...
                                                                                                  continuation_options.time_step.options);
                [beta_k_tmp, ~, this.last_mean_time_grid, this.last_mean_time_step_diagnostics, mean_prior_diagnostics] = ...
                    pt_cont.Pseudo_Time_Continuation_Forward_Euler_Backtracking(theta_traj, [], ...
                                                                                backtracking_options);
                this.last_mean_time_step_prior_diagnostics = mean_prior_diagnostics;
            else
                theta_traj = MD_Discrepancy_Parameter_Trajectory(continuation_options.num_steps, sample_idx);
                beta_k_tmp = pt_cont.Pseudo_Time_Continuation_Forward_Euler(theta_traj);
                if sample_idx == 0
                    this.last_mean_time_step_prior_diagnostics = struct();
                    this.last_mean_time_grid = linspace(0, 1, continuation_options.num_steps + 1);
                    this.last_mean_time_step_diagnostics = {};
                end
            end

            beta_k = beta_k_tmp(:, end);
            u_k = pt_cont.sen_op.current_u;
            z_k = pt_cont.sen_op.current_z;

        end

        function [u_k, z_k, beta_k, prior_diagnostics, time_grid, step_diagnostics] = Posterior_Update_Sample_From_Mean_Core(this, beta_bar, sample_idx, sen_op)
            if nargin < 4 || isempty(sen_op)
                sen_op = MD_Continuation_Sensitivity_Operators(this.post_sampling, this.hessian_analysis, this.options.samples.discard_cache, this.options.samples.materialize_breve_samples);
            end

            qn_prec = MD_Quasi_Newton_Preconditioner(this.hessian_analysis);
            pt_cont = Pseudo_Time_Continuation(beta_bar, sen_op, qn_prec);
            prior_diagnostics = struct();
            time_grid = [];
            step_diagnostics = {};
            sample_options = this.options.samples;

            if sample_options.time_step.use
                theta_traj = MD_Discrepancy_Mean_To_Sample_Trajectory(1, sample_idx);

                backtracking_options = sample_options.time_step.options;
                backtracking_options.time_step_callback = ...
                    @(beta, theta, time_index, beta_dot) sen_op.Sample_A_Priori_Time_Step_Heuristic(beta, theta, time_index, beta_dot, sample_idx, ...
                                                                                                    sample_options.time_step.epsilon_nl, ...
                                                                                                    sample_options.time_step.options);
                [beta_k_tmp, ~, time_grid, step_diagnostics, prior_diagnostics] = ...
                    pt_cont.Pseudo_Time_Continuation_Forward_Euler_Backtracking(theta_traj, [], ...
                                                                                backtracking_options);
            else
                theta_traj = MD_Discrepancy_Mean_To_Sample_Trajectory(sample_options.num_mean_to_sample_steps, sample_idx);
                beta_k_tmp = pt_cont.Pseudo_Time_Continuation_Forward_Euler(theta_traj);
                time_grid = linspace(0, 1, sample_options.num_mean_to_sample_steps + 1);
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
            sen_op = MD_Continuation_Sensitivity_Operators(this.post_sampling, this.hessian_analysis, this.options.samples.discard_cache, this.options.samples.materialize_breve_samples);
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

            sen_op = MD_Continuation_Sensitivity_Operators(this.post_sampling, this.hessian_analysis);

            u = this.opt_prob_interface.State_Solve(z);
            delta = sen_op.Discrepancy_Evaluation_Mean(z);

            [val, grad_u, grad_z] = this.opt_prob_interface.Objective_Function(u + delta, z);
            z_tmp1 = sen_op.Apply_Discrepancy_z_Jacobian_Transpose_Mean(grad_u);
            z_tmp2 = this.opt_prob_interface.Apply_Solution_Operator_z_Jacobian_Transpose(grad_u, z);
            grad = grad_z + z_tmp1 + z_tmp2;

        end

    end

end
