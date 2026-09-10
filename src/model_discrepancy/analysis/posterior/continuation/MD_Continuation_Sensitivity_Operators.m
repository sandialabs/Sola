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

        % Adaptive posterior-sample continuation data.
        %
        % breve_R satisfies
        %
        %   breve_R * breve_R' approx Sigma_beta.
        %
        % Each entry of breve_samplers stores one persistent lazy matrix-normal
        % realization for the corresponding posterior sample index.
        breve_R
        breve_samplers
        lazy_sampling_tol
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

            % Correction for nonlinear state map S(z).
            %
            % Apply_RS_Hessian uses grad_u evaluated at S(z).  The continuation
            % objective uses grad_u evaluated at S(z) + t*d(z).  Therefore we need
            %
            %   S_zz(z)[z_in]^* * (grad_u_corrected - grad_u_low_fidelity).
            %
            grad_u_corrected = this.opt_prob_interface.Misfit_Gradient(this.current_u + delta, this.current_z);
            grad_u_low_fidelity = this.opt_prob_interface.Misfit_Gradient(this.current_u, this.current_z);
            grad_u_diff = grad_u_corrected - grad_u_low_fidelity;
            z_out = z_out + this.opt_prob_interface.Apply_Solution_Operator_z_Hessian_Adjoint(z_in, grad_u_diff, this.current_z);

            beta_out = this.hessian_analysis.Apply_V_Transpose(z_out);

            u_D = this.current_disc_ops.Apply_Beta_Jacobian(beta_in, this.current_t);
            x = this.opt_prob_interface.Apply_Misfit_Hessian(u_D, this.current_u + delta, this.current_z);
            z_S_adj = this.opt_prob_interface.Apply_Solution_Operator_z_Jacobian_Transpose(x, this.current_z);
            beta_out = beta_out + this.hessian_analysis.Apply_V_Transpose(z_S_adj);

            beta_out = beta_out + this.current_disc_ops.Apply_Beta_Jacobian_Transpose(x, this.current_t);

            u_S = this.opt_prob_interface.Apply_Solution_Operator_z_Jacobian(z_in, this.current_z);
            x = this.opt_prob_interface.Apply_Misfit_Hessian(u_S, this.current_u + delta, this.current_z);
            beta_out = beta_out + this.current_disc_ops.Apply_Beta_Jacobian_Transpose(x, this.current_t);
        end

        function [beta_out] = Apply_B(this, beta, theta_traj, time_index)

            beta = beta(:);

            this.State_Evaluation(beta, theta_traj, time_index);

            delta = this.current_disc_ops.Eval(beta, this.current_t);

            D = this.current_disc_ops.Apply_Theta_Jacobian(beta);
            x = this.opt_prob_interface.Apply_Misfit_Hessian(D, this.current_u + delta, this.current_z);
            z_S_adj = this.opt_prob_interface.Apply_Solution_Operator_z_Jacobian_Transpose(x, this.current_z);
            beta_out = this.hessian_analysis.Apply_V_Transpose(z_S_adj);
            beta_out = beta_out + this.current_disc_ops.Apply_Beta_Jacobian_Transpose(x, this.current_t);
            state_grad = this.opt_prob_interface.Misfit_Gradient(this.current_u + delta, this.current_z);
            beta_out = beta_out + this.current_disc_ops.Apply_Beta_Theta_Hessian(state_grad);

        end

    end

    methods

        function this = MD_Continuation_Sensitivity_Operators(post_sampling, hessian_analysis)
            arguments
                post_sampling MD_Posterior_Sampling
                hessian_analysis MD_Hessian_Analysis
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

            this.Mz_Wz_inv_Mz_Z_minus_z_opt = this.post_data.Mz_Wz_inv_Mz_Z - this.post_data.Mz_Wz_inv_Mz_z_opt;
            this.Mz_Wz_inv_Mz_yi = 0 * this.Mz_Wz_inv_Mz_Z_minus_z_opt;
            this.si = zeros(1, this.post_data.N);

            for i = 1:this.post_data.N
                this.Mz_Wz_inv_Mz_yi(:, i) = this.post_data.Mz_Wz_inv_Mz_Z * this.post_data.g_vecs(:, i) - ...
                    sum(this.post_data.g_vecs(:, i)) * this.post_data.Mz_Wz_inv_Mz_z_opt;
                this.si(i) = sum(this.post_data.g_vecs(:, i)) - this.z_opt' * this.Mz_Wz_inv_Mz_yi(:, i);
            end

            this.lazy_sampling_tol = 1e-10;
            this.breve_R = this.Compute_Breve_Beta_Covariance_Factor();
            this.breve_samplers = cell(this.post_data.num_samples, 1);
        end

        function [] = State_Evaluation(this, beta, theta_traj, time_index)
            beta = beta(:);

            t = theta_traj.Get_Time(time_index);
            sample_idx = theta_traj.Get_Sample_Index();

            state_changed = max([ ...
                abs(t - this.current_t), ...
                norm(beta - this.current_beta), ...
                abs(sample_idx - this.current_sample_idx)]) > 1.e-15;

            if state_changed

                this.current_t = t;
                this.current_beta = beta;
                this.current_sample_idx = sample_idx;
                this.current_z = this.z_opt + this.hessian_analysis.Apply_V(beta);
                this.current_u = this.opt_prob_interface.State_Solve(this.current_z);
                this.current_disc_ops = this.Get_Continuation_Beta_Discrepancy_Ops(sample_idx);

            end

            if sample_idx > 0
                this.Get_Breve_Sampler(sample_idx);
            end

        end

        function [disc_ops] = Get_Continuation_Beta_Discrepancy_Ops(this, sample_idx)
            % Return beta-space discrepancy operators for continuation.
            %
            % All returned operators are consistent with
            %
            %   z_beta = z_opt + V beta
            %
            % and with the continuation state
            %
            %   S(z_beta) + t D(beta).
            %
            % The methods are:
            %
            %   Eval(beta, t)
            %       Returns t * D(beta).
            %
            %   Apply_Beta_Jacobian(beta_in, t)
            %       Returns t * D_beta beta_in.
            %
            %   Apply_Beta_Jacobian_Transpose(u_in, t)
            %       Returns t * D_beta' u_in.
            %
            %   Apply_Theta_Jacobian(beta)
            %       Returns D(beta), i.e. derivative of t*D(beta) with respect to t.
            %
            %   Apply_Beta_Theta_Hessian(u_in)
            %       Returns D_beta' u_in, i.e. derivative with respect to beta
            %       of Apply_Theta_Jacobian paired against u_in.

            assert(sample_idx >= 0 && sample_idx <= this.post_data.num_samples && floor(sample_idx) == sample_idx, ...
                'sample_idx must be an integer in [0, num_samples].');

            if sample_idx == 0
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

        % ------------------------------------------------------------
        % Beta-space breve covariance and sampler cache
        % ------------------------------------------------------------

        function [R, Sigma_beta] = Compute_Breve_Beta_Covariance_Factor(this)

            if isempty(this.hessian_analysis.evals)
                r = length(this.z_opt);
            else
                r = length(this.hessian_analysis.evals);
            end

            V = zeros(length(this.z_opt), r);

            for j = 1:r
                e = zeros(r, 1);
                e(j) = 1.0;
                V(:, j) = this.hessian_analysis.Apply_V(e);
            end

            Mz_V = this.z_prior_interface.Apply_M_z(V);
            Wz_inv_Mz_V = this.z_prior_interface.Apply_W_z_Inverse(Mz_V);

            if isempty(this.post_data.Zc_Mz_Wz_inv_Mz_Zc)
                tmp_rhs = Wz_inv_Mz_V;
            else
                tmp_rhs = Wz_inv_Mz_V - this.post_data.Wz_inv_Mz_Zc * linsolve( ...
                        this.post_data.Zc_Mz_Wz_inv_Mz_Zc, this.post_data.Mz_Zc' * Wz_inv_Mz_V);
            end

            Sigma_beta = Mz_V' * tmp_rhs;
            Sigma_beta = 0.5 * (Sigma_beta + Sigma_beta');

            [U, D] = eig(Sigma_beta);
            lambda = real(diag(D));

            lambda_scale = max(abs(lambda));
            if isempty(lambda_scale)
                lambda_scale = 0;
            end

            tol = 1e-12 * max(1.0, lambda_scale);

            keep = lambda > tol;

            if any(lambda < -tol)
                disp('Warning: Compute_Breve_Beta_Covariance_Factor found negative eigenvalues below tolerance.');
            end

            if any(keep)
                R = U(:, keep) * diag(sqrt(lambda(keep)));
            else
                R = zeros(r, 0);
            end

        end

        function sampler = Get_Breve_Sampler(this, sample_idx)

            assert(sample_idx >= 1 && sample_idx <= this.post_data.num_samples && floor(sample_idx) == sample_idx, ...
                   'sample_idx must be an integer in [1, num_samples] for posterior samples.');

            if isempty(this.breve_samplers{sample_idx})
                this.breve_samplers{sample_idx} = MD_Breve_Beta_Sampler(this.breve_R, this.post_sampling.u_prior_interface, length(this.u_opt), this.lazy_sampling_tol);
            end

            sampler = this.breve_samplers{sample_idx};

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
            u_out = (1 / this.post_data.alpha_d) * u_out;
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
        % Posterior sample explicit kernels excluding breve term
        % ------------------------------------------------------------

        function [u_out] = Discrepancy_Evaluation_Sample_Explicit(this, z, sample_idx)

            assert(sample_idx >= 1 && sample_idx <= this.post_data.num_samples && floor(sample_idx) == sample_idx, ...
                   'sample_idx must be an integer in [1, num_samples].');

            u_out = this.Discrepancy_Evaluation_Mean(z);
            dz = z - this.z_opt;

            u_hat = zeros(size(u_out));

            for i = 1:this.post_data.N
                sgi = sum(this.post_data.g_vecs(:, i));
                coeff = (1 / sqrt(this.post_data.Mu(i, i))) * (sgi + this.Mz_Wz_inv_Mz_yi(:, i)' * dz);
                u_hat = u_hat + coeff * this.post_data.ui_hat{i}(:, sample_idx);
            end

            u_hat = sqrt(this.post_data.alpha_d) * u_hat;
            u_out = u_out + u_hat;

        end

        function [u_out] = Apply_Discrepancy_z_Jacobian_Sample_Explicit(this, z_in, sample_idx)

            assert(sample_idx >= 1 && sample_idx <= this.post_data.num_samples && floor(sample_idx) == sample_idx, ...
                   'sample_idx must be an integer in [1, num_samples].');

            u_out = this.Apply_Discrepancy_z_Jacobian_Mean(z_in);
            u_hat = zeros(size(u_out));

            for i = 1:this.post_data.N
                coeff = (1 / sqrt(this.post_data.Mu(i, i))) * (this.Mz_Wz_inv_Mz_yi(:, i)' * z_in);
                u_hat = u_hat + coeff * this.post_data.ui_hat{i}(:, sample_idx);
            end

            u_hat = sqrt(this.post_data.alpha_d) * u_hat;
            u_out = u_out + u_hat;

        end

        function [z_out] = Apply_Discrepancy_z_Jacobian_Transpose_Sample_Explicit(this, u_in, sample_idx)

            assert(sample_idx >= 1 && sample_idx <= this.post_data.num_samples && floor(sample_idx) == sample_idx, ...
                   'sample_idx must be an integer in [1, num_samples].');

            z_out = this.Apply_Discrepancy_z_Jacobian_Transpose_Mean(u_in);
            z_hat = zeros(size(z_out));

            for i = 1:this.post_data.N
                ui_hat_idx = this.post_data.ui_hat{i}(:, sample_idx);
                coeff = (1 / sqrt(this.post_data.Mu(i, i))) * (ui_hat_idx' * u_in);
                z_hat = z_hat + coeff * this.Mz_Wz_inv_Mz_yi(:, i);
            end

            z_hat = sqrt(this.post_data.alpha_d) * z_hat;
            z_out = z_out + z_hat;

        end

        % ------------------------------------------------------------
        % Posterior sample beta-space explicit-plus-breve kernels
        % ------------------------------------------------------------

        function [u_out] = Discrepancy_Evaluation_Sample_Beta(this, beta, sample_idx)

            assert(sample_idx >= 1 && sample_idx <= this.post_data.num_samples && floor(sample_idx) == sample_idx, ...
                   'sample_idx must be an integer in [1, num_samples].');

            beta = beta(:);
            z = this.z_opt + this.hessian_analysis.Apply_V(beta);
            u_out = this.Discrepancy_Evaluation_Sample_Explicit(z, sample_idx);
            sampler = this.Get_Breve_Sampler(sample_idx);
            u_out = u_out + sampler.Eval(beta);

        end

        function [u_out] = Apply_Discrepancy_Beta_Jacobian_Sample(this, beta_in, sample_idx)

            assert(sample_idx >= 1 && sample_idx <= this.post_data.num_samples && floor(sample_idx) == sample_idx, ...
                   'sample_idx must be an integer in [1, num_samples].');

            beta_in = beta_in(:);
            z_in = this.hessian_analysis.Apply_V(beta_in);
            u_out = this.Apply_Discrepancy_z_Jacobian_Sample_Explicit(z_in, sample_idx);
            sampler = this.Get_Breve_Sampler(sample_idx);
            u_out = u_out + sampler.Apply_Jacobian(beta_in);

        end

        function [beta_out] = Apply_Discrepancy_Beta_Jacobian_Transpose_Sample(this, u_in, sample_idx)

            assert(sample_idx >= 1 && sample_idx <= this.post_data.num_samples && floor(sample_idx) == sample_idx, ...
                   'sample_idx must be an integer in [1, num_samples].');

            z_out = this.Apply_Discrepancy_z_Jacobian_Transpose_Sample_Explicit(u_in, sample_idx);
            beta_out = this.hessian_analysis.Apply_V_Transpose(z_out);
            sampler = this.Get_Breve_Sampler(sample_idx);
            beta_out = beta_out + sampler.Apply_Jacobian_Transpose(u_in);

        end

    end

end