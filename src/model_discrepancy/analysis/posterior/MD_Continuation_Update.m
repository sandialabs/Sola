classdef MD_Continuation_Update < handle

    properties
        md_post_sampling
        md_hessian_analysis
        opt_prob_interface
        u_opt
        z_opt
        num_continuation_steps
        step_size

        Mz_Wz_inv_Mz_Z_minus_z_opt
        Mz_Wz_inv_Mz_yi
        si
    end

    methods

        function this = MD_Continuation_Update(md_post_sampling, md_hessian_analysis, num_continuation_steps)
            this.md_post_sampling = md_post_sampling;
            this.md_hessian_analysis = md_hessian_analysis;
            this.opt_prob_interface = md_hessian_analysis.opt_prob_interface;
            this.u_opt = md_post_sampling.data_interface.u_opt;
            this.z_opt = md_post_sampling.data_interface.z_opt;

            this.Mz_Wz_inv_Mz_Z_minus_z_opt = this.md_post_sampling.post_data.Mz_Wz_inv_Mz_Z - this.md_post_sampling.post_data.Mz_Wz_inv_Mz_z_opt;
            this.Mz_Wz_inv_Mz_yi = 0 * this.Mz_Wz_inv_Mz_Z_minus_z_opt;
            for i = 1:this.md_post_sampling.post_data.N
                this.Mz_Wz_inv_Mz_yi(:, i) = this.md_post_sampling.post_data.Mz_Wz_inv_Mz_Z * this.md_post_sampling.post_data.g_vecs(:, i) - sum(this.md_post_sampling.post_data.g_vecs(:, i)) * this.md_post_sampling.post_data.Mz_Wz_inv_Mz_z_opt;
                this.si(i) = sum(this.md_post_sampling.post_data.g_vecs(:, i)) - this.z_opt' * this.Mz_Wz_inv_Mz_yi(:, i);
            end

            this.num_continuation_steps = num_continuation_steps;
            this.step_size = 1 / num_continuation_steps;
        end

        function [u, z, beta] = Posterior_Update_Mean(this)
            u = zeros(length(this.u_opt), this.num_continuation_steps + 1);
            z = zeros(length(this.z_opt), this.num_continuation_steps + 1);
            t = linspace(0, 1, this.num_continuation_steps + 1);
            beta = zeros(length(this.md_hessian_analysis.evals), this.num_continuation_steps + 1);

            u(:, 1) = this.u_opt;
            z(:, 1) = this.z_opt;

            for k = 1:this.num_continuation_steps
                % Predictive step for beta
                Btheta_n = this.Apply_B_beta(u(:, k), beta(:, k), t(k));
                beta_pert = -this.Apply_Parameterized_RS_Hessian_Inverse_beta(Btheta_n, u(:, k), beta(:, k), t(k));
                beta_pred = beta(:, k) + this.step_size * beta_pert;

                % Predictive step for z and u
                z_pred = this.z_opt + this.md_hessian_analysis.evecs * beta_pred;
                u_pred = this.opt_prob_interface.State_Solve(z_pred);

                % Corrective step for beta
                Jbeta_val = this.md_hessian_analysis.evecs' * this.Gradient_J_z(u_pred, z_pred, t(k + 1));
                beta(:, k + 1) = beta_pred - this.Apply_Parameterized_RS_Hessian_Inverse_beta(Jbeta_val, u_pred, beta_pred, t(k + 1));

                % Corrective step for z and u
                z(:, k + 1) = this.z_opt + this.md_hessian_analysis.evecs * beta(:, k + 1);
                u(:, k + 1) = this.opt_prob_interface.State_Solve(z(:, k + 1));
            end
        end

        function [val, grad] = Jhat_Posterior_beta(this, beta)
            z = this.z_opt + this.md_hessian_analysis.evecs * beta;
            [val, J_grad_z] = this.Jhat_Posterior(z);
            grad = this.md_hessian_analysis.evecs' * J_grad_z;
        end

        function [val, grad] = Jhat_Posterior(this, z)
            u = this.opt_prob_interface.State_Solve(z);
            delta = this.Discrepancy_Evaluation(z, 1);

            [val, grad_u, grad_z] = this.opt_prob_interface.Objective_Function(u + delta, z);
            z_tmp1 = this.Apply_Discrepancy_z_Jacobian_transpose(grad_u, 1);
            z_tmp2 = this.opt_prob_interface.Apply_Solution_Operator_z_Jacobian_Transpose(grad_u, z);
            grad = grad_z + z_tmp1 + z_tmp2;
        end

        function [Btheta_n] = Apply_B_beta(this, u_n, beta_n, t_n)
            z_n = this.z_opt + this.md_hessian_analysis.evecs * beta_n;
            Btheta_n = this.md_hessian_analysis.evecs' * this.Apply_B(u_n, z_n, t_n);
        end

        function [beta_out] = Apply_Parameterized_RS_Hessian_Inverse_beta(this, beta_in, u_n, beta_n, t_n)
            beta_out = 0 * beta_in;
            for k = 1:size(beta_in, 2)
                tol = 1.e-7;
                max_iter = length(beta_n) + 10;
                [beta_out(:, k), flag, relres, iter, resvec] = pcg(@(x)this.Apply_Parameterized_RS_Hessian_beta(x, u_n, beta_n, t_n), beta_in(:, k), tol, max_iter);
                if flag ~= 0
                    disp(['CG did not converge; flag: ', num2str(flag), ', relres: ', num2str(relres)]);
                end
            end
        end

        function [beta_out] = Apply_Parameterized_RS_Hessian_beta(this, beta_in, u_n, beta_n, t_n)
            z_n = this.z_opt + this.md_hessian_analysis.evecs * beta_n;
            beta_out = this.md_hessian_analysis.evecs' * this.Apply_Parameterized_RS_Hessian(this.md_hessian_analysis.evecs * beta_in, u_n, z_n, t_n);
        end

        function [z_out] = Gradient_J_z(this, u_n, z_n, t_n)
            delta = this.Discrepancy_Evaluation(z_n, t_n);
            [~, grad_u, grad_z] = this.opt_prob_interface.Objective_Function(u_n + delta, z_n);
            z_tmp1 = this.Apply_Discrepancy_z_Jacobian_transpose(grad_u, t_n);
            z_tmp2 = this.opt_prob_interface.Apply_Solution_Operator_z_Jacobian_Transpose(grad_u, z_n);
            z_out = grad_z + z_tmp1 + z_tmp2;
        end

        function [Btheta_n] = Apply_B(this, u_n, z_n, t_n)
            delta = this.Discrepancy_Evaluation(z_n, t_n);

            u_tmp1 = this.Apply_Discrepancy_theta_Jacobian(z_n);
            u_tmp2 = this.opt_prob_interface.Apply_Misfit_Hessian(u_tmp1, u_n + delta, z_n);
            z_tmp1 = this.opt_prob_interface.Apply_Solution_Operator_z_Jacobian_Transpose(u_tmp2, z_n);

            z_tmp2 = this.Apply_Discrepancy_z_Jacobian_transpose(u_tmp2, t_n);

            state_grad = this.opt_prob_interface.Misfit_Gradient(u_n + delta, z_n);
            z_tmp3 = this.Apply_Discrepancy_z_theta_Mean(state_grad);

            Btheta_n = z_tmp1 + z_tmp2 + z_tmp3;
        end

        function [z_out] = Apply_Parameterized_RS_Hessian(this, z_in, u_n, z_n, t_n)
            delta = this.Discrepancy_Evaluation(z_n, t_n);

            z_out = this.opt_prob_interface.Apply_RS_Hessian(z_in, z_n); % NOTE: this computes J_{zz} + S_z' * J_uu * S_z
            z_out = z_out(:);

            u_tmp1 = this.Apply_Discrepancy_z_Jacobian(z_in, t_n);
            u_tmp2 = this.opt_prob_interface.Apply_Misfit_Hessian(u_tmp1, u_n + delta, z_n);
            z_out = z_out + this.opt_prob_interface.Apply_Solution_Operator_z_Jacobian_Transpose(u_tmp2, z_n);

            z_out = z_out + this.Apply_Discrepancy_z_Jacobian_transpose(u_tmp2, t_n);

            u_tmp3 = this.opt_prob_interface.Apply_Solution_Operator_z_Jacobian(z_in, z_n);
            u_tmp4 = this.opt_prob_interface.Apply_Misfit_Hessian(u_tmp3, u_n + delta, z_n);
            z_out = z_out + this.Apply_Discrepancy_z_Jacobian_transpose(u_tmp4, t_n);
        end

        function [u_out] = Discrepancy_Evaluation(this, z_n, t_n)
            N = this.md_post_sampling.post_data.N;
            u_out = 0 * this.u_opt;
            for ell = 1:N
                coeff = this.md_post_sampling.post_data.a_ell(ell) + z_n' * this.Mz_Wz_inv_Mz_Z_minus_z_opt(:, ell);
                u_out = u_out + coeff * this.md_post_sampling.post_data.u_ell(:, ell);
                for i = 1:N
                    coeff = this.md_post_sampling.post_data.b_i_ell(i, ell) * (this.si(i) + z_n' * this.Mz_Wz_inv_Mz_yi(:, i));
                    u_out = u_out - coeff * this.md_post_sampling.post_data.u_i_ell{i}(:, ell);
                end
            end
            u_out = (t_n / this.md_post_sampling.post_data.alpha_d) * u_out;
        end

        function [u_out] = Apply_Discrepancy_z_Jacobian(this, z_in, t_n)
            N = this.md_post_sampling.post_data.N;
            u = 0 * this.u_opt;
            for ell = 1:N
                coeff = this.Mz_Wz_inv_Mz_Z_minus_z_opt(:, ell)' * z_in;
                u = u + coeff * this.md_post_sampling.post_data.u_ell(:, ell);
                for i = 1:N
                    coeff = this.md_post_sampling.post_data.b_i_ell(i, ell) * (this.Mz_Wz_inv_Mz_yi(:, i)' * z_in);
                    u = u - coeff * this.md_post_sampling.post_data.u_i_ell{i}(:, ell);
                end
            end

            u_out = t_n * (1 / this.md_post_sampling.post_data.alpha_d) * u;
        end

        function [z_out] = Apply_Discrepancy_z_Jacobian_transpose(this, u_in, t_n)
            N = this.md_post_sampling.post_data.N;
            z = 0 * this.z_opt;
            for ell = 1:N
                z = z + (this.md_post_sampling.post_data.u_ell(:, ell)' * u_in) * this.Mz_Wz_inv_Mz_Z_minus_z_opt(:, ell);
                for i = 1:N
                    coeff = this.md_post_sampling.post_data.b_i_ell(i, ell) * (this.md_post_sampling.post_data.u_i_ell{i}(:, ell)' * u_in);
                    vec = this.Mz_Wz_inv_Mz_yi(:, i);
                    z = z - coeff * vec;
                end
            end

            z_out = t_n * (1 / this.md_post_sampling.post_data.alpha_d) * z;
        end

        function [u_out] = Apply_Discrepancy_theta_Jacobian(this, z_n)
            % Simplifies since discrepancy is linear in theta
            u_out = this.Discrepancy_Evaluation(z_n, 1.0);
        end

        function [z_out] = Apply_Discrepancy_z_theta_Mean(this, u_in)
            % Simplifies since discrepancy is linear in theta
            z_out = this.Apply_Discrepancy_z_Jacobian_transpose(u_in, 1.0);
        end

        % -------------------------------------------
        % IN PROGRESS: Continuation for samples
        % -------------------------------------------

        % Steps that remain:
        % * Implement the two remaining procedures: Apply_Discrepancy_z_Jacobian_transpose_Sample & Apply_Discrepancy_z_Jacobian_Sample
        % * After proper testing, refactor existing code so that sample_idx = 0 corresponds to mean.
        % * Clean up naming convention and documentation for clarity

        function [u_ks, z_ks, beta_ks] = Posterior_Update_Samples(this)
            num_samples = this.md_post_sampling.post_data.num_samples;
            u_ks = zeros(length(this.u_opt), num_samples);
            z_ks = zeros(length(this.z_opt), num_samples);
            beta_ks = zeros(length(this.md_hessian_analysis.eval), num_samples);
            for sample_idx = 1:num_samples
                [u_k, z_k, beta_k] = this.Posterior_Update_Sample(sample_idx);
                u_ks(:, sample_idx) = u_k;
                z_ks(:, sample_idx) = z_k;
                beta_ks(:, sample_idx) = beta_k;
            end
        end

        % [WIP]
        function [u_k, z_k, beta_k] = Posterior_Update_Sample(this, sample_idx)
            t = linspace(0, 1, this.num_continuation_steps + 1);
            u_k = this.u_opt;
            beta_k = zeros(length(this.md_hessian_analysis.evals));

            for k = 1:this.num_continuation_steps
                % Predictive step for beta
                Btheta_n = this.Apply_B_beta_Sample(u_k, beta_k, t(k), sample_idx);
                beta_pert = -this.Apply_Parameterized_RS_Hessian_Inverse_beta_Sample(Btheta_n, u_k, beta_k, t(k), sample_idx);
                beta_pred = beta_k + this.step_size * beta_pert;

                % Predictive step for z and u
                z_pred = this.z_opt + this.md_hessian_analysis.evecs * beta_pred;
                u_pred = this.opt_prob_interface.State_Solve(z_pred);

                % Corrective step for beta
                Jbeta_val = this.md_hessian_analysis.evecs' * this.Gradient_J_z_Sample(u_pred, z_pred, t(k + 1), sample_idx);
                beta_new = beta_pred - this.Apply_Parameterized_RS_Hessian_Inverse_beta_Sample(Jbeta_val, u_pred, beta_pred, t(k + 1), sample_idx);

                % Update indices
                beta_k = beta_new;
                z_k = this.z_opt + this.md_hessian_analysis.evecs * beta_new;
                u_k = this.opt_prob_interface.State_Solve(z_k);
            end
        end

        % [Verify]
        function [val, grad] = Jhat_Posterior_beta_Sample(this, beta, sample_idx)
            z = this.z_opt + this.md_hessian_analysis.evecs * beta;
            [val, J_grad_z] = this.Jhat_Posterior_Sample(z, sample_idx);
            grad = this.md_hessian_analysis.evecs' * J_grad_z;
        end

        % [Verify]
        function [val, grad] = Jhat_Posterior_Sample(this, z, sample_idx)
            u = this.opt_prob_interface.State_Solve(z);
            delta = this.Discrepancy_Evaluation_Sample(z, 1, sample_idx);

            [val, grad_u, grad_z] = this.opt_prob_interface.Objective_Function(u + delta, z);
            z_tmp1 = this.Apply_Discrepancy_z_Jacobian_transpose_Sample(grad_u, 1, sample_idx);
            z_tmp2 = this.opt_prob_interface.Apply_Solution_Operator_z_Jacobian_Transpose(grad_u, z);
            grad = grad_z + z_tmp1 + z_tmp2;
        end

        % [Verify]
        function [Btheta_n] = Apply_B_beta_Sample(this, u_n, beta_n, t_n, sample_idx)
            z_n = this.z_opt + this.md_hessian_analysis.evecs * beta_n;
            Btheta_n = this.md_hessian_analysis.evecs' * this.Apply_B_Sample(u_n, z_n, t_n, sample_idx);
        end

        % [Verify]
        function [beta_out] = Apply_Parameterized_RS_Hessian_Inverse_beta_Sample(this, beta_in, u_n, beta_n, t_n, sample_idx)
            beta_out = 0 * beta_in;
            for k = 1:size(beta_in, 2)
                tol = 1.e-7;
                max_iter = length(beta_n) + 10;
                [beta_out(:, k), flag, relres, iter, resvec] = pcg(@(x)this.Apply_Parameterized_RS_Hessian_beta_Sample(x, u_n, beta_n, t_n, sample_idx), beta_in(:, k), tol, max_iter);
                if flag ~= 0
                    disp(['CG did not converge; flag: ', num2str(flag), ', relres: ', num2str(relres)]);
                end
            end
        end

        % [Verify]
        function [beta_out] = Apply_Parameterized_RS_Hessian_beta_Sample(this, beta_in, u_n, beta_n, t_n, sample_idx)
            z_n = this.z_opt + this.md_hessian_analysis.evecs * beta_n;
            beta_out = this.md_hessian_analysis.evecs' * this.Apply_Parameterized_RS_Hessian_Sample(this.md_hessian_analysis.evecs * beta_in, u_n, z_n, t_n, sample_idx);
        end

        % [Verify]
        function [z_out] = Gradient_J_z_Sample(this, u_n, z_n, t_n, sample_idx)
            delta = this.Discrepancy_Evaluation_Sample(z_n, t_n, sample_idx);
            [~, grad_u, grad_z] = this.opt_prob_interface.Objective_Function_Sample(u_n + delta, z_n, sample_idx);
            z_tmp1 = this.Apply_Discrepancy_z_Jacobian_transpose_Sample(grad_u, t_n, sample_idx);
            z_tmp2 = this.opt_prob_interface.Apply_Solution_Operator_z_Jacobian_Transpose(grad_u, z_n);
            z_out = grad_z + z_tmp1 + z_tmp2;
        end

        % [Verify]
        function [Btheta_n] = Apply_B_Sample(this, u_n, z_n, t_n, sample_idx)
            delta = this.Discrepancy_Evaluation_Sample(z_n, t_n, sample_idx);

            u_tmp1 = this.Apply_Discrepancy_theta_Jacobian_Sample(z_n, sample_idx);
            u_tmp2 = this.opt_prob_interface.Apply_Misfit_Hessian(u_tmp1, u_n + delta, z_n);
            z_tmp1 = this.opt_prob_interface.Apply_Solution_Operator_z_Jacobian_Transpose(u_tmp2, z_n);

            z_tmp2 = this.Apply_Discrepancy_z_Jacobian_transpose_Sample(u_tmp2, t_n, sample_idx);

            state_grad = this.opt_prob_interface.Misfit_Gradient(u_n + delta, z_n);
            z_tmp3 = this.Apply_Discrepancy_z_theta_Sample(state_grad, sample_idx);

            Btheta_n = z_tmp1 + z_tmp2 + z_tmp3;
        end

        % [Verify]
        function [z_out] = Apply_Parameterized_RS_Hessian_Sample(this, z_in, u_n, z_n, t_n, sample_idx)
            delta = this.Discrepancy_Evaluation_Sample(z_n, t_n, sample_idx);

            z_out = this.opt_prob_interface.Apply_RS_Hessian(z_in, z_n); % NOTE: this computes J_{zz} + S_z' * J_uu * S_z
            z_out = z_out(:);

            u_tmp1 = this.Apply_Discrepancy_z_Jacobian_Sample(z_in, t_n, sample_idx);
            u_tmp2 = this.opt_prob_interface.Apply_Misfit_Hessian(u_tmp1, u_n + delta, z_n);
            z_out = z_out + this.opt_prob_interface.Apply_Solution_Operator_z_Jacobian_Transpose(u_tmp2, z_n);

            z_out = z_out + this.Apply_Discrepancy_z_Jacobian_transpose_Sample(u_tmp2, t_n, sample_idx);

            u_tmp3 = this.opt_prob_interface.Apply_Solution_Operator_z_Jacobian(z_in, z_n);
            u_tmp4 = this.opt_prob_interface.Apply_Misfit_Hessian(u_tmp3, u_n + delta, z_n);
            z_out = z_out + this.Apply_Discrepancy_z_Jacobian_transpose_Sample(u_tmp4, t_n, sample_idx);
        end

        % TODO: Rework this to avoid computing everything
        function [u_out] = Discrepancy_Evaluation_Sample(this, z_n, t_n, sample_idx)
            [~, delta_samples] = this.md_post_sampling.Posterior_Discrepancy_Sample(z_n);
            u_out = t_n * delta_samples{1};
            u_out = u_out(:, sample_idx);
        end

        % TODO: IMPLEMENT THIS ROUTINE
        function [u_out] = Apply_Discrepancy_z_Jacobian_Sample(this, z_in, t_n, sample_idx)
            u_out_mean = this.Apply_Discrepancy_z_Jacobian(z_in, t_n);
            u = 0 * this.u_opt; % TODO: Implement
            u_out = u_out_mean + t_n * u;
        end

        % TODO: IMPLEMENT THIS ROUTINE
        function [z_out] = Apply_Discrepancy_z_Jacobian_transpose_Sample(this, u_in, t_n, sample_idx)
            z_out_mean = this.Apply_Discrepancy_z_Jacobian_transpose(u_in, t_n);
            z = 0 * this.z_opt; % TODO: Implement
            z_out = z_out_mean + t_n * z;
        end

        function [u_out] = Apply_Discrepancy_theta_Jacobian_Sample(this, z_n, sample_idx)
            % Simplifies since discrepancy is linear in theta
            u_out = this.Discrepancy_Evaluation_Sample(z_n, 1.0, sample_idx);
        end

        function [z_out] = Apply_Discrepancy_z_theta_Sample(this, u_in)
            % Simplifies since discrepancy is linear in theta
            z_out = this.Apply_Discrepancy_z_Jacobian_transpose_Sample(u_in, 1.0, sample_idx);
        end

    end

end
