%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef MD_Continuation_Sensitivity_Operators < Sensitivity_Operators

    properties
        post_sampling
        post_data
        hessian_analysis
        opt_prob_interface
        data_interface
        z_prior_interface
        u_opt
        z_opt

        Mz_Wz_inv_Mz_Z_minus_z_opt
        Mz_Wz_inv_Mz_yi
        si

        current_t
        current_u
        current_beta
        current_z
        current_disc_ops
        current_sample_idx
        current_is_mean_to_sample

        % Posterior-sample continuation data.
        %
        % Each entry of breve_samplers stores one persistent breve
        % matrix-normal realization for the corresponding posterior sample
        % index.  By default, the realization is matrix-free/lazy.  If
        % materialize_breve_samples is true, each breve sample is revealed on
        % all beta basis vectors when it is first constructed, yielding a
        % fixed dense output_dim-by-r beta-space matrix for deterministic
        % trial evaluations.
        breve_samplers
        discard_cache
        lazy_sampling_tol
        materialize_breve_samples
    end

    methods (Access = public)

        function [grad, val] = Gradient(this, beta, theta_traj, time_index)

            beta = beta(:);

            this.State_Evaluation(beta, theta_traj, time_index);

            delta = this.current_disc_ops.Eval(beta, this.current_t);

            [val, grad_u, grad_z] = this.opt_prob_interface.Objective_Function(this.current_u + delta, this.current_z);
            beta_grad_z = this.hessian_analysis.Apply_V_Transpose(grad_z);
            z_S_adj = this.opt_prob_interface.Apply_Solution_Operator_z_Jacobian_Transpose(grad_u, this.current_z);
            beta_grad_S = this.hessian_analysis.Apply_V_Transpose(z_S_adj);
            beta_grad_D = this.current_disc_ops.Apply_Beta_Jacobian_Transpose(grad_u, this.current_t);
            grad = beta_grad_z + beta_grad_S + beta_grad_D;

        end

        function [beta_out] = Apply_Hessian(this, beta_in, beta, theta_traj, time_index)
            beta_in = beta_in(:);
            beta = beta(:);

            this.State_Evaluation(beta, theta_traj, time_index);
            delta = this.current_disc_ops.Eval(beta, this.current_t);
            z_in = this.hessian_analysis.Apply_V(beta_in);

            % Base reduced-space Hessian contribution.
            z_out = this.opt_prob_interface.Apply_RS_Hessian(z_in, this.current_z);

            % Corrections to the low-fidelity reduced Hessian.
            %
            % Apply_RS_Hessian differentiates J(S(z), z) instead of J(S(z) + delta(z,t), z).
            % Therefore, the low-fidelity reduced Hessian must be corrected for
            % both changes in grad_u and J_uu.
            % (1) S_zz(z)[z_in]^* * (grad_u_corrected - grad_u_lofi)
            % (2) S_z^*(z)[(Juu u_S)_corrected - (Juu u_S)_lofi]
            u_S = this.opt_prob_interface.Apply_Solution_Operator_z_Jacobian(z_in, this.current_z);
            grad_u_corrected = this.opt_prob_interface.Misfit_Gradient(this.current_u + delta, this.current_z);
            grad_u_lofi = this.opt_prob_interface.Misfit_Gradient(this.current_u, this.current_z);
            grad_u_diff = grad_u_corrected - grad_u_lofi;
            z_out = z_out + this.opt_prob_interface.Apply_Solution_Operator_z_Hessian_Adjoint(z_in, grad_u_diff, this.current_z);

            x_corrected = this.opt_prob_interface.Apply_Misfit_Hessian(u_S, this.current_u + delta, this.current_z);
            x_lofi = this.opt_prob_interface.Apply_Misfit_Hessian(u_S, this.current_u, this.current_z);
            z_S_adj = this.opt_prob_interface.Apply_Solution_Operator_z_Jacobian_Transpose(x_corrected - x_lofi, this.current_z);
            z_out = z_out + z_S_adj;

            beta_out = this.hessian_analysis.Apply_V_Transpose(z_out);

            u_D = this.current_disc_ops.Apply_Beta_Jacobian(beta_in, this.current_t);
            x = this.opt_prob_interface.Apply_Misfit_Hessian(u_D, this.current_u + delta, this.current_z);
            z_S_adj = this.opt_prob_interface.Apply_Solution_Operator_z_Jacobian_Transpose(x, this.current_z);
            beta_out = beta_out + this.hessian_analysis.Apply_V_Transpose(z_S_adj);

            beta_out = beta_out + this.current_disc_ops.Apply_Beta_Jacobian_Transpose(x, this.current_t);

            x = this.opt_prob_interface.Apply_Misfit_Hessian(u_S, this.current_u + delta, this.current_z);
            beta_out = beta_out + this.current_disc_ops.Apply_Beta_Jacobian_Transpose(x, this.current_t);
        end

        function [beta_out] = Apply_B(this, beta, theta_traj, time_index)

            beta = beta(:);

            this.State_Evaluation(beta, theta_traj, time_index);

            delta = this.current_disc_ops.Eval(beta, this.current_t);

            u_tmp = this.current_disc_ops.Apply_Theta_Jacobian(beta);
            u_tmp = this.opt_prob_interface.Apply_Misfit_Hessian(u_tmp, this.current_u + delta, this.current_z);
            z_tmp = this.opt_prob_interface.Apply_Solution_Operator_z_Jacobian_Transpose(u_tmp, this.current_z);
            beta_out = this.hessian_analysis.Apply_V_Transpose(z_tmp);
            beta_out = beta_out + this.current_disc_ops.Apply_Beta_Jacobian_Transpose(u_tmp, this.current_t);
            state_grad = this.opt_prob_interface.Misfit_Gradient(this.current_u + delta, this.current_z);
            beta_out = beta_out + this.current_disc_ops.Apply_Beta_Theta_Hessian(state_grad);

        end


    end

    methods

        function this = MD_Continuation_Sensitivity_Operators(post_sampling, hessian_analysis, discard_cache, materialize_breve_samples)
            arguments
                post_sampling MD_Posterior_Sampling
                hessian_analysis MD_Hessian_Analysis
                discard_cache (1,1) logical = true
                materialize_breve_samples (1,1) logical = false
            end

            this.post_sampling = post_sampling;
            this.post_data = post_sampling.post_data;
            this.hessian_analysis = hessian_analysis;
            this.opt_prob_interface = hessian_analysis.opt_prob_interface;
            this.data_interface = post_sampling.data_interface;
            this.z_prior_interface = post_sampling.z_prior_interface;
            this.u_opt = post_sampling.data_interface.u_opt;
            this.z_opt = post_sampling.data_interface.z_opt;

            this.current_t = inf;
            this.current_beta = inf;
            this.current_sample_idx = inf;
            this.current_is_mean_to_sample = false;

            this.Mz_Wz_inv_Mz_Z_minus_z_opt = this.post_data.Mz_Wz_inv_Mz_Z - this.post_data.Mz_Wz_inv_Mz_z_opt;
            this.Mz_Wz_inv_Mz_yi = 0 * this.Mz_Wz_inv_Mz_Z_minus_z_opt;
            this.si = zeros(1, this.post_data.N);

            for i = 1:this.post_data.N
                this.Mz_Wz_inv_Mz_yi(:, i) = this.post_data.Mz_Wz_inv_Mz_Z * this.post_data.g_vecs(:, i) - ...
                    sum(this.post_data.g_vecs(:, i)) * this.post_data.Mz_Wz_inv_Mz_z_opt;
                this.si(i) = sum(this.post_data.g_vecs(:, i)) - this.z_opt' * this.Mz_Wz_inv_Mz_yi(:, i);
            end

            this.lazy_sampling_tol = 1e-10;
            this.discard_cache = discard_cache;
            this.materialize_breve_samples = materialize_breve_samples;
            this.breve_samplers = cell(this.post_data.num_samples, 1);
        end

        function [dt, diagnostics] = Mean_A_Priori_Time_Step_Heuristic(this, beta, theta_traj, time_index, beta_dot, epsilon_nl, options)
            arguments
                this
                beta (:, 1) {mustBeNumeric}
                theta_traj Auxillary_Parameter_Trajectory
                time_index (1, 1) {mustBeNumeric}
                beta_dot (:, 1) {mustBeNumeric}
                epsilon_nl (1, 1) {mustBeNumeric, mustBeNonnegative} = 0.1
                options struct = struct()
            end

            apply_increment_jacobian = @(beta_in) this.Apply_Discrepancy_z_Jacobian_Mean(this.hessian_analysis.Apply_V(beta_in));
            apply_increment_jacobian_transpose = @(u_in) this.hessian_analysis.Apply_V_Transpose(this.Apply_Discrepancy_z_Jacobian_Transpose_Mean(u_in));

            [dt, diagnostics] = this.A_Priori_Time_Step_Heuristic(beta, theta_traj, time_index, beta_dot, ...
                                                                  apply_increment_jacobian, ...
                                                                  apply_increment_jacobian_transpose, ...
                                                                  epsilon_nl, options);
        end

        function [dt, diagnostics] = Sample_A_Priori_Time_Step_Heuristic(this, beta, theta_traj, time_index, beta_dot, sample_idx, epsilon_nl, options)
            arguments
                this
                beta (:, 1) {mustBeNumeric}
                theta_traj Auxillary_Parameter_Trajectory
                time_index (1, 1) {mustBeNumeric}
                beta_dot (:, 1) {mustBeNumeric}
                sample_idx (1, 1) {mustBeNumeric}
                epsilon_nl (1, 1) {mustBeNumeric, mustBeNonnegative} = 0.1
                options struct = struct()
            end

            apply_increment_jacobian = @(beta_in) this.Apply_Discrepancy_Beta_Jacobian_Sample_Minus_Mean(beta_in, sample_idx);
            apply_increment_jacobian_transpose = @(u_in) this.Apply_Discrepancy_Beta_Jacobian_Transpose_Sample_Minus_Mean(u_in, sample_idx);

            [dt, diagnostics] = this.A_Priori_Time_Step_Heuristic(beta, theta_traj, time_index, beta_dot, ...
                                                                  apply_increment_jacobian, ...
                                                                  apply_increment_jacobian_transpose, ...
                                                                  epsilon_nl, options);
        end

        function [dt, diagnostics] = A_Priori_Time_Step_Heuristic(this, beta, theta_traj, time_index, beta_dot, apply_increment_jacobian, apply_increment_jacobian_transpose, epsilon_nl, options)
            arguments
                this
                beta (:, 1) {mustBeNumeric}
                theta_traj Auxillary_Parameter_Trajectory
                time_index (1, 1) {mustBeNumeric}
                beta_dot (:, 1) {mustBeNumeric}
                apply_increment_jacobian function_handle
                apply_increment_jacobian_transpose function_handle
                epsilon_nl (1, 1) {mustBeNumeric, mustBeNonnegative} = 0.1
                options struct = struct()
            end

            this.State_Evaluation(beta, theta_traj, time_index);

            fd_rel_step = this.Get_Option(options, 'fd_rel_step', 1.e-4);
            min_denominator = this.Get_Option(options, 'min_denominator', 1.e-14);
            min_dt = this.Get_Option(options, 'min_dt', 1.e-3);
            max_dt = this.Get_Option(options, 'max_dt', 1.0);
            epsilon_curv = this.Get_Option(options, 'epsilon_curv', 1.0);

            z_dot = this.hessian_analysis.Apply_V(beta_dot);
            z_dot_norm = norm(z_dot);

            if z_dot_norm <= min_denominator
                A = 0.0;
                first_norm = 0.0;
                second_norm = 0.0;
                fd_step = 0.0;
            else
                fd_step = fd_rel_step * max(1.0, norm(this.current_z)) / z_dot_norm;
                u_plus = this.opt_prob_interface.State_Solve(this.current_z + fd_step * z_dot);
                u_minus = this.opt_prob_interface.State_Solve(this.current_z - fd_step * z_dot);
                second_order_action = (u_plus - 2 * this.current_u + u_minus) / (fd_step^2);
                first_order_action = this.opt_prob_interface.Apply_Solution_Operator_z_Jacobian(z_dot, this.current_z);

                second_norm = this.Weighted_u_Norm(second_order_action);
                first_norm = this.Weighted_u_Norm(first_order_action);

                if first_norm <= min_denominator
                    A = 0.0;
                else
                    A = second_norm / (2 * first_norm);
                end
            end

            eta = this.Relative_Gauss_Newton_Curvature(beta, theta_traj, time_index, ...
                                                       apply_increment_jacobian, ...
                                                       apply_increment_jacobian_transpose, ...
                                                       min_denominator);

            if A <= min_denominator
                dt_nl = inf;
            else
                dt_nl = epsilon_nl / A;
            end

            if eta <= min_denominator
                dt_curv = inf;
            else
                dt_curv = sqrt(epsilon_curv / eta);
            end

            dt = min([max_dt, 1.0, dt_nl, dt_curv]);
            if ~isfinite(dt) || dt <= 0.0
                dt = min(max_dt, 1.0);
            end
            dt = min(max_dt, max(min_dt, dt));

            diagnostics = struct();
            diagnostics.t = this.current_t;
            diagnostics.A = A;
            diagnostics.eta = eta;
            diagnostics.dt_nonlinearity = dt_nl;
            diagnostics.dt_curvature = dt_curv;
            diagnostics.epsilon_nl = epsilon_nl;
            diagnostics.epsilon_curv = epsilon_curv;
            diagnostics.first_order_norm = first_norm;
            diagnostics.second_order_norm = second_norm;
            diagnostics.fd_step = fd_step;
            diagnostics.dt = dt;
        end

        function [eta] = Relative_Gauss_Newton_Curvature(this, beta, theta_traj, time_index, apply_increment_jacobian, apply_increment_jacobian_transpose, regularization)
            beta = beta(:);
            r = length(beta);
            H_ref = zeros(r, r);
            H_delta = zeros(r, r);

            delta = this.current_disc_ops.Eval(beta, this.current_t);
            for j = 1:r
                e_j = zeros(r, 1);
                e_j(j) = 1.0;
                H_ref(:, j) = this.Apply_Hessian(e_j, beta, theta_traj, time_index);

                u_j = apply_increment_jacobian(e_j);
                R_u_j = this.opt_prob_interface.Apply_Misfit_Hessian(u_j, this.current_u + delta, this.current_z);
                H_delta(:, j) = apply_increment_jacobian_transpose(R_u_j);
            end

            H_ref = 0.5 * (H_ref + H_ref');
            H_delta = 0.5 * (H_delta + H_delta');
            H_ref = H_ref + regularization * eye(r);

            if norm(H_delta, 'fro') <= regularization
                eta = 0.0;
                return;
            end

            try
                vals = eig(H_delta, H_ref, 'vector');
            catch
                vals = eig(pinv(H_ref) * H_delta, 'vector');
            end
            vals = real(vals(isfinite(vals)));
            if isempty(vals)
                eta = 0.0;
            else
                eta = max(0.0, max(vals));
            end
        end

        function [nrm] = Weighted_u_Norm(this, u)
            try
                Wu_u = this.post_sampling.u_prior_interface.Apply_W_u(u);
                if isempty(Wu_u)
                    nrm = norm(u);
                else
                    nrm = sqrt(abs(u' * Wu_u));
                end
            catch
                nrm = norm(u);
            end
        end

        function [value] = Get_Option(this, options, name, default_value) %#ok<INUSL>
            if isfield(options, name) && ~isempty(options.(name))
                value = options.(name);
            else
                value = default_value;
            end
        end

        function [] = State_Evaluation(this, beta, theta_traj, time_index)
            beta = beta(:);

            t = theta_traj.Get_Time(time_index);
            sample_idx = theta_traj.Get_Sample_Index();
            is_mean_to_sample = this.Is_Mean_To_Sample_Trajectory(theta_traj);

            if max([abs(t - this.current_t), norm(beta - this.current_beta), abs(sample_idx - this.current_sample_idx)]) > 1.e-15 || ...
                    is_mean_to_sample ~= this.current_is_mean_to_sample

                this.current_t = t;
                this.current_beta = beta;
                this.current_sample_idx = sample_idx;
                this.current_is_mean_to_sample = is_mean_to_sample;
                this.current_z = this.z_opt + this.hessian_analysis.Apply_V(beta);
                this.current_u = this.opt_prob_interface.State_Solve(this.current_z);
                if is_mean_to_sample
                    this.current_disc_ops = this.Get_Mean_To_Sample_Beta_Discrepancy_Ops(sample_idx);
                else
                    this.current_disc_ops = this.Get_Continuation_Beta_Discrepancy_Ops(sample_idx);
                end

            end

            if sample_idx > 0
                this.Get_Breve_Sampler(sample_idx);
            end

        end

        function [is_mean_to_sample] = Is_Mean_To_Sample_Trajectory(this, theta_traj)
            if ismethod(theta_traj, 'Is_Mean_To_Sample')
                is_mean_to_sample = theta_traj.Is_Mean_To_Sample();
            else
                is_mean_to_sample = false;
            end
        end

        function [disc_ops] = Get_Continuation_Beta_Discrepancy_Ops(this, sample_idx)
            assert(sample_idx >= 0 && sample_idx <= this.post_data.num_samples && floor(sample_idx) == sample_idx, ...
                'sample_idx must be an integer in [0, num_samples].');

            if sample_idx == 0 % Mean
                disc_ops.Eval = @(beta, t) t * this.Discrepancy_Evaluation_Mean(this.z_opt + this.hessian_analysis.Apply_V(beta));
                disc_ops.Apply_Beta_Jacobian = @(beta_in, t) t * this.Apply_Discrepancy_z_Jacobian_Mean(this.hessian_analysis.Apply_V(beta_in));
                disc_ops.Apply_Beta_Jacobian_Transpose = @(u_in, t) t * this.hessian_analysis.Apply_V_Transpose(this.Apply_Discrepancy_z_Jacobian_Transpose_Mean(u_in));
                disc_ops.Apply_Theta_Jacobian = @(beta) this.Discrepancy_Evaluation_Mean(this.z_opt + this.hessian_analysis.Apply_V(beta));
                disc_ops.Apply_Beta_Theta_Hessian = @(u_in) this.hessian_analysis.Apply_V_Transpose(this.Apply_Discrepancy_z_Jacobian_Transpose_Mean(u_in));

            else
                this.Get_Breve_Sampler(sample_idx);
                disc_ops.Eval = @(beta, t) t * this.Discrepancy_Evaluation_Sample_Beta(beta, sample_idx);
                disc_ops.Apply_Beta_Jacobian = @(beta_in, t) t * this.Apply_Discrepancy_Beta_Jacobian_Sample(beta_in, sample_idx);
                disc_ops.Apply_Beta_Jacobian_Transpose = @(u_in, t) t * this.Apply_Discrepancy_Beta_Jacobian_Transpose_Sample(u_in, sample_idx);
                disc_ops.Apply_Theta_Jacobian = @(beta) this.Discrepancy_Evaluation_Sample_Beta(beta, sample_idx);
                disc_ops.Apply_Beta_Theta_Hessian = @(u_in) this.Apply_Discrepancy_Beta_Jacobian_Transpose_Sample(u_in, sample_idx);
            end
        end

        function [disc_ops] = Get_Mean_To_Sample_Beta_Discrepancy_Ops(this, sample_idx)
            assert(sample_idx >= 1 && sample_idx <= this.post_data.num_samples && floor(sample_idx) == sample_idx, ...
                'sample_idx must be an integer in [1, num_samples].');

            this.Get_Breve_Sampler(sample_idx);
            disc_ops.Eval = @(beta, s) this.Discrepancy_Evaluation_Mean_To_Sample_Beta(beta, s, sample_idx);
            disc_ops.Apply_Beta_Jacobian = @(beta_in, s) this.Apply_Discrepancy_Beta_Jacobian_Mean_To_Sample(beta_in, s, sample_idx);
            disc_ops.Apply_Beta_Jacobian_Transpose = @(u_in, s) this.Apply_Discrepancy_Beta_Jacobian_Transpose_Mean_To_Sample(u_in, s, sample_idx);
            disc_ops.Apply_Theta_Jacobian = @(beta) this.Discrepancy_Evaluation_Sample_Minus_Mean_Beta(beta, sample_idx);
            disc_ops.Apply_Beta_Theta_Hessian = @(u_in) this.Apply_Discrepancy_Beta_Jacobian_Transpose_Sample_Minus_Mean(u_in, sample_idx);
        end

        % ------------------------------------------------------------
        % Beta-space breve covariance and sampler cache
        % ------------------------------------------------------------

        function sampler = Get_Breve_Sampler(this, sample_idx)

            if isempty(this.breve_samplers{sample_idx})
                if this.discard_cache
                    this.Clear_Breve_Samplers_Except(sample_idx);
                end

                this.breve_samplers{sample_idx} = MD_Breve_Beta_Sampler(this.hessian_analysis, ...
                    this.z_prior_interface, this.post_data, this.post_sampling.u_prior_interface, ...
                    length(this.u_opt), this.lazy_sampling_tol, this.materialize_breve_samples);
            end

            sampler = this.breve_samplers{sample_idx};

        end

        function [] = Clear_Breve_Samplers_Except(this, sample_idx)

            for i = 1:length(this.breve_samplers)
                if i ~= sample_idx
                    this.breve_samplers{i} = [];
                end
            end

        end

        % ------------------------------------------------------------
        % Discrepancy kernels: Mean
        % ------------------------------------------------------------

        function [u_out] = Discrepancy_Evaluation_Mean(this, z)
            N = this.post_data.N;
            u_out = zeros(size(this.data_interface.u_opt));
            for ell = 1:N
                coeff = this.post_data.a_ell(ell) + z' * this.Mz_Wz_inv_Mz_Z_minus_z_opt(:, ell);
                u_out = u_out + coeff * this.post_data.u_ell(:, ell);
                for i = 1:N
                    coeff = this.post_data.b_i_ell(i, ell) * (this.si(i) + z' * this.Mz_Wz_inv_Mz_yi(:, i));
                    u_out = u_out - coeff * this.post_data.u_i_ell{i}(:, ell);
                end
            end
            u_out = (1 / this.post_data.alpha_d) * u_out + this.data_interface.data_shift;
        end

        function [u_out] = Apply_Discrepancy_z_Jacobian_Mean(this, z_in)
            N = this.post_data.N;
            u = zeros(size(this.data_interface.u_opt));
            for ell = 1:N
                u = u + (this.Mz_Wz_inv_Mz_Z_minus_z_opt(:, ell)' * z_in) * this.post_data.u_ell(:, ell);
                for i = 1:N
                    coeff = this.post_data.b_i_ell(i, ell) * (this.Mz_Wz_inv_Mz_yi(:, i)' * z_in);
                    u = u - coeff * this.post_data.u_i_ell{i}(:, ell);
                end
            end

            u_out = (1 / this.post_data.alpha_d) * u;
        end

        function [z_out] = Apply_Discrepancy_z_Jacobian_Transpose_Mean(this, u_in)
            N = this.post_data.N;
            z = zeros(size(this.z_opt));
            for ell = 1:N
                z = z + (this.post_data.u_ell(:, ell)' * u_in) * this.Mz_Wz_inv_Mz_Z_minus_z_opt(:, ell);
                for i = 1:N
                    coeff = this.post_data.b_i_ell(i, ell) * (this.post_data.u_i_ell{i}(:, ell)' * u_in);
                    z = z - coeff * this.Mz_Wz_inv_Mz_yi(:, i);
                end
            end

            z_out = (1 / this.post_data.alpha_d) * z;
        end

        % ------------------------------------------------------------
        % Posterior sample kernels
        % ------------------------------------------------------------

        function [u_out] = Discrepancy_Evaluation_Mean_To_Sample_Beta(this, beta, s, sample_idx)

            u_mean = this.Discrepancy_Evaluation_Mean(this.z_opt + this.hessian_analysis.Apply_V(beta));
            u_diff = this.Discrepancy_Evaluation_Sample_Beta(beta, sample_idx) - u_mean;
            u_out = u_mean + s * u_diff;

        end

        function [u_out] = Apply_Discrepancy_Beta_Jacobian_Mean_To_Sample(this, beta_in, s, sample_idx)

            u_mean = this.Apply_Discrepancy_z_Jacobian_Mean(this.hessian_analysis.Apply_V(beta_in));
            u_diff = this.Apply_Discrepancy_Beta_Jacobian_Sample(beta_in, sample_idx) - u_mean;
            u_out = u_mean + s * u_diff;

        end

        function [beta_out] = Apply_Discrepancy_Beta_Jacobian_Transpose_Mean_To_Sample(this, u_in, s, sample_idx)

            beta_mean = this.hessian_analysis.Apply_V_Transpose(this.Apply_Discrepancy_z_Jacobian_Transpose_Mean(u_in));
            beta_diff = this.Apply_Discrepancy_Beta_Jacobian_Transpose_Sample(u_in, sample_idx) - beta_mean;
            beta_out = beta_mean + s * beta_diff;

        end

        function [u_out] = Discrepancy_Evaluation_Sample_Minus_Mean_Beta(this, beta, sample_idx)

            u_out = this.Discrepancy_Evaluation_Sample_Beta(beta, sample_idx) - ...
                this.Discrepancy_Evaluation_Mean(this.z_opt + this.hessian_analysis.Apply_V(beta));

        end

        function [u_out] = Apply_Discrepancy_Beta_Jacobian_Sample_Minus_Mean(this, beta_in, sample_idx)

            u_out = this.Apply_Discrepancy_Beta_Jacobian_Sample(beta_in, sample_idx) - ...
                this.Apply_Discrepancy_z_Jacobian_Mean(this.hessian_analysis.Apply_V(beta_in));

        end

        function [beta_out] = Apply_Discrepancy_Beta_Jacobian_Transpose_Sample_Minus_Mean(this, u_in, sample_idx)

            beta_out = this.Apply_Discrepancy_Beta_Jacobian_Transpose_Sample(u_in, sample_idx) - ...
                this.hessian_analysis.Apply_V_Transpose(this.Apply_Discrepancy_z_Jacobian_Transpose_Mean(u_in));

        end

        function [u_out] = Discrepancy_Evaluation_Sample_Beta(this, beta, sample_idx)

            beta = beta(:);
            dz = this.hessian_analysis.Apply_V(beta);
            
            % Mean component
            u_out = this.Discrepancy_Evaluation_Mean(this.z_opt + dz);

            % Hat component
            u_hat = zeros(size(u_out));
            for i = 1:this.post_data.N
                sgi = sum(this.post_data.g_vecs(:, i));
                coeff = (1 / sqrt(this.post_data.Mu(i, i))) * (sgi + this.Mz_Wz_inv_Mz_yi(:, i)' * dz);
                u_hat = u_hat + coeff * this.post_data.ui_hat{i}(:, sample_idx);
            end
            u_out = u_out + sqrt(this.post_data.alpha_d) * u_hat;
            
            % Breve component
            sampler = this.Get_Breve_Sampler(sample_idx);
            u_out = u_out + sampler.Eval(beta);

        end

        function [u_out] = Apply_Discrepancy_Beta_Jacobian_Sample(this, beta_in, sample_idx)

            beta_in = beta_in(:);
            z_in = this.hessian_analysis.Apply_V(beta_in);

            % Mean component
            u_out = this.Apply_Discrepancy_z_Jacobian_Mean(z_in);

            % Hat component
            u_hat = zeros(size(u_out));
            for i = 1:this.post_data.N
                coeff = (1 / sqrt(this.post_data.Mu(i, i))) * (this.Mz_Wz_inv_Mz_yi(:, i)' * z_in);
                u_hat = u_hat + coeff * this.post_data.ui_hat{i}(:, sample_idx);
            end
            u_out = u_out + sqrt(this.post_data.alpha_d) * u_hat;
            
            % Breve component
            sampler = this.Get_Breve_Sampler(sample_idx);
            u_out = u_out + sampler.Apply_Jacobian(beta_in);

        end

        function [beta_out] = Apply_Discrepancy_Beta_Jacobian_Transpose_Sample(this, u_in, sample_idx)

            % Mean component
            z_out = this.Apply_Discrepancy_z_Jacobian_Transpose_Mean(u_in);
            
            % Hat component
            z_hat = zeros(size(z_out));
            for i = 1:this.post_data.N
                ui_hat_idx = this.post_data.ui_hat{i}(:, sample_idx);
                coeff = (1 / sqrt(this.post_data.Mu(i, i))) * (ui_hat_idx' * u_in);
                z_hat = z_hat + coeff * this.Mz_Wz_inv_Mz_yi(:, i);
            end
            z_out = z_out + sqrt(this.post_data.alpha_d) * z_hat;
            beta_out = this.hessian_analysis.Apply_V_Transpose(z_out);
            
            % Breve component
            sampler = this.Get_Breve_Sampler(sample_idx);
            beta_out = beta_out + sampler.Apply_Jacobian_Transpose(u_in);

        end

    end

end
