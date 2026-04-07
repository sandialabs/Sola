classdef MD_Continuation_Update < handle

    properties
        md_post_sampling
        md_hessian_analysis
        opt_prob_interface
        u_opt
        z_opt
        num_continuation_steps

        Mz_Wz_inv_Mz_Z_minus_z_opt
        Mz_Wz_inv_Mz_yi
        si
        r % NEW
    end

    methods

        function this = MD_Continuation_Update(md_post_sampling, md_hessian_analysis, num_continuation_steps)
            this.md_post_sampling = md_post_sampling;
            this.md_hessian_analysis = md_hessian_analysis;
            this.opt_prob_interface = md_hessian_analysis.opt_prob_interface;
            this.u_opt = md_post_sampling.data_interface.u_opt;
            this.z_opt = md_post_sampling.data_interface.z_opt;
            if isempty(md_hessian_analysis.evals)
                this.r = length(this.z_opt);
            else
                this.r = length(md_hessian_analysis.evals);
            end
            this.num_continuation_steps = num_continuation_steps;
        end

        % ------------------------------------------------------------
        % High-level functions
        % ------------------------------------------------------------

        function [u, z, beta] = Posterior_Update_Mean(this)
            discOps = this.Get_Discrepancy_Ops(0);
            [u, z, beta] = this.Posterior_Update_Core(discOps);
        end

        function [u_ks, z_ks, beta_ks] = Posterior_Update_Samples(this)
            num_samples = this.md_post_sampling.post_data.num_samples;
            u_ks = zeros(length(this.u_opt), num_samples);
            z_ks = zeros(length(this.z_opt), num_samples);
            beta_ks = zeros(this.r, num_samples);

            for sample_idx = 1:num_samples
                discOps = this.Get_Discrepancy_Ops(sample_idx);
                [u_k, z_k, beta_k] = this.Posterior_Update_Core(discOps);
                u_ks(:, sample_idx) = u_k;
                z_ks(:, sample_idx) = z_k;
                beta_ks(:, sample_idx) = beta_k;
            end
        end

        % ------------------------------------------------------------
        % Core continuation driver (shared by mean and samples)
        % ------------------------------------------------------------

        function [u_k, z_k, beta_k] = Posterior_Update_Core(this, discOps)
            t = linspace(0, 1, this.num_continuation_steps + 1);

            u_k = this.u_opt;
            beta_k = zeros(this.r, 1);
            z_k = this.z_opt;

            for k = 1:this.num_continuation_steps
                % Predictive step for beta
                Btheta_n = this.Apply_B_beta(u_k, beta_k, t(k), discOps);
                beta_pert = -this.Apply_Parameterized_RS_Hessian_Inverse_beta(Btheta_n, u_k, beta_k, t(k), discOps);
                beta_pred = beta_k + (1 / this.num_continuation_steps) * beta_pert;

                % Predictive step for z and u
                z_pred = this.z_opt + this.md_hessian_analysis.V_apply(beta_pred);
                u_pred = this.opt_prob_interface.State_Solve(z_pred);

                % Corrective step for beta
                Jbeta_val = this.md_hessian_analysis.V_transpose_apply(this.Gradient_J_z(u_pred, z_pred, t(k + 1), discOps));
                beta_new = beta_pred - this.Apply_Parameterized_RS_Hessian_Inverse_beta(Jbeta_val, u_pred, beta_pred, t(k + 1), discOps);

                % Update state
                beta_k = beta_new;
                z_k = this.z_opt + this.md_hessian_analysis.V_apply(beta_new);
                u_k = this.opt_prob_interface.State_Solve(z_k);
            end
        end

        % ------------------------------------------------------------
        % Objective helpers (optional; kept for completeness)
        % ------------------------------------------------------------

        function [val, grad] = Jhat_Posterior_beta(this, beta, discOps)
            z = this.z_opt + this.md_hessian_analysis.V_apply(beta);
            [val, J_grad_z] = this.Jhat_Posterior(z, discOps);
            grad = this.md_hessian_analysis.V_transpose_apply(J_grad_z);
        end

        function [val, grad] = Jhat_Posterior(this, z, discOps)
            u = this.opt_prob_interface.State_Solve(z);
            delta = discOps.Eval(z, 1);

            [val, grad_u, grad_z] = this.opt_prob_interface.Objective_Function(u + delta, z);
            z_tmp1 = discOps.ApplyJzT(z, grad_u, 1);
            z_tmp2 = this.opt_prob_interface.Apply_Solution_Operator_z_Jacobian_Transpose(grad_u, z);
            grad = grad_z + z_tmp1 + z_tmp2;
        end

        % ------------------------------------------------------------
        % Continuation building blocks (shared)
        % ------------------------------------------------------------

        function [Btheta_n] = Apply_B_beta(this, u_n, beta_n, t_n, discOps)
            z_n = this.z_opt + this.md_hessian_analysis.V_apply(beta_n);
            Btheta_n = this.md_hessian_analysis.V_transpose_apply(this.Apply_B(u_n, z_n, t_n, discOps));
        end

        function [beta_out] = Apply_Parameterized_RS_Hessian_Inverse_beta(this, beta_in, u_n, beta_n, t_n, discOps)
            beta_out = 0 * beta_in;
            for k = 1:size(beta_in, 2)
                tol = 1.e-7;
                max_iter = 2 * length(beta_n) + 10;
                [beta_out(:, k), flag, relres] = pcg( ...
                                                     @(x)this.Apply_Parameterized_RS_Hessian_beta(x, u_n, beta_n, t_n, discOps), ...
                                                     beta_in(:, k), tol, max_iter);

                if flag ~= 0
                    disp(['CG did not converge; flag: ', num2str(flag), ', relres: ', num2str(relres)]);
                end
            end
        end

        function [beta_out] = Apply_Parameterized_RS_Hessian_beta(this, beta_in, u_n, beta_n, t_n, discOps)
            z_n = this.z_opt + this.md_hessian_analysis.V_apply(beta_n);
            beta_out = this.md_hessian_analysis.V_transpose_apply(this.Apply_Parameterized_RS_Hessian(this.md_hessian_analysis.V_apply(beta_in), u_n, z_n, t_n, discOps));
        end

        function [z_out] = Gradient_J_z(this, u_n, z_n, t_n, discOps)
            delta = discOps.Eval(z_n, t_n);
            [~, grad_u, grad_z] = this.opt_prob_interface.Objective_Function(u_n + delta, z_n);
            z_tmp1 = discOps.ApplyJzT(z_n, grad_u, t_n);
            z_tmp2 = this.opt_prob_interface.Apply_Solution_Operator_z_Jacobian_Transpose(grad_u, z_n);
            z_out = grad_z + z_tmp1 + z_tmp2;
        end

        function [Btheta_n] = Apply_B(this, u_n, z_n, t_n, discOps)
            delta = discOps.Eval(z_n, t_n);

            u_tmp1 = discOps.ApplyJtheta(z_n);
            u_tmp2 = this.opt_prob_interface.Apply_Misfit_Hessian(u_tmp1, u_n + delta, z_n);
            z_tmp1 = this.opt_prob_interface.Apply_Solution_Operator_z_Jacobian_Transpose(u_tmp2, z_n);

            z_tmp2 = discOps.ApplyJzT(z_n, u_tmp2, t_n);

            state_grad = this.opt_prob_interface.Misfit_Gradient(u_n + delta, z_n);
            z_tmp3 = discOps.ApplyJzthetaT(z_n, state_grad);

            Btheta_n = z_tmp1 + z_tmp2 + z_tmp3;
        end

        function [z_out] = Apply_Parameterized_RS_Hessian(this, z_in, u_n, z_n, t_n, discOps)
            delta = discOps.Eval(z_n, t_n);

            % NOTE: computes J_{zz} + S_z' * J_uu * S_z
            z_out = this.opt_prob_interface.Apply_RS_Hessian(z_in, z_n);
            z_out = z_out(:);

            % Term involving delta_z(z_n) * z_in
            u_tmp1 = discOps.ApplyJz(z_n, z_in, t_n);
            u_tmp2 = this.opt_prob_interface.Apply_Misfit_Hessian(u_tmp1, u_n + delta, z_n);
            z_out = z_out + this.opt_prob_interface.Apply_Solution_Operator_z_Jacobian_Transpose(u_tmp2, z_n);

            z_out = z_out + discOps.ApplyJzT(z_n, u_tmp2, t_n);

            % Cross term: S_z * z_in then delta_z^T
            u_tmp3 = this.opt_prob_interface.Apply_Solution_Operator_z_Jacobian(z_in, z_n);
            u_tmp4 = this.opt_prob_interface.Apply_Misfit_Hessian(u_tmp3, u_n + delta, z_n);
            z_out = z_out + discOps.ApplyJzT(z_n, u_tmp4, t_n);
        end

        % ------------------------------------------------------------
        % Discrepancy ops factory
        % ------------------------------------------------------------

        function discOps = Get_Discrepancy_Ops(this, sample_idx)
            num_samples = this.md_post_sampling.post_data.num_samples;
            assert(sample_idx >= 0 && sample_idx <= num_samples && floor(sample_idx) == sample_idx, ...
                   'sample_idx must be an integer in [0, num_samples].');

            if sample_idx == 0
                discOps.Eval      = @(z, t) t * this.md_post_sampling.Discrepancy_Evaluation_Mean(z);
                discOps.ApplyJz   = @(z, z_in, t) t * this.md_post_sampling.Apply_Discrepancy_z_Jacobian_Mean(z_in);
                discOps.ApplyJzT  = @(z, u,  t) t * this.md_post_sampling.Apply_Discrepancy_z_Jacobian_transpose_Mean(u);

                % Linear-in-theta simplifications (mean)
                discOps.ApplyJtheta   = @(z) this.md_post_sampling.Discrepancy_Evaluation_Mean(z);
                discOps.ApplyJzthetaT = @(z, u) this.md_post_sampling.Apply_Discrepancy_z_Jacobian_transpose_Mean(u);
            else
                discOps.Eval      = @(z, t) t * this.md_post_sampling.Discrepancy_Evaluation_Sample(z, sample_idx);
                discOps.ApplyJz   = @(z, z_in, t) t * this.md_post_sampling.Apply_Discrepancy_z_Jacobian_Sample(z, z_in, sample_idx);
                discOps.ApplyJzT  = @(z, u,  t) t * this.md_post_sampling.Apply_Discrepancy_z_Jacobian_transpose_Sample(z, u, sample_idx);

                % Linear-in-theta simplifications (sample)
                discOps.ApplyJtheta   = @(z) this.md_post_sampling.Discrepancy_Evaluation_Sample(z, sample_idx);
                discOps.ApplyJzthetaT = @(z, u) this.md_post_sampling.Apply_Discrepancy_z_Jacobian_transpose_Sample(z, u, sample_idx);
            end
        end

    end

end
