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
            % Useful for performing finite-difference checks
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
            u_out = (1 / this.md_post_sampling.post_data.alpha_d) * u_out;
        end

        function [z_out] = Apply_Discrepancy_z_theta_Mean(this, u_in)
            z_out = this.Apply_Discrepancy_z_Jacobian_transpose(u_in, 1.0);
        end

        % REMOVED FUNCTIONS (Reference for translation purposes)

        % function [u, z] = Posterior_Update_Mean(this)
        %     u = zeros(length(this.u_opt), this.num_continuation_steps + 1);
        %     z = zeros(length(this.z_opt), this.num_continuation_steps + 1);
        %     t = linspace(0, 1, this.num_continuation_steps + 1);

        %     u(:, 1) = this.u_opt;
        %     z(:, 1) = this.z_opt;

        %     for k = 1:this.num_continuation_steps
        %         Btheta_n = this.Apply_B(u(:, k), z(:, k), t(k));
        %         z_pert = -this.Apply_Parameterized_RS_Hessian_Inverse(Btheta_n, u(:, k), z(:, k), t(k));
        %         z(:, k + 1) = z(:, k) + this.step_size * z_pert;
        %         u(:, k + 1) = this.opt_prob_interface.State_Solve(z(:, k + 1));
        %     end
        % end

        % function [u, z] = Posterior_Update_Mean_PC(this)
        %     u = zeros(length(this.u_opt), this.num_continuation_steps + 1);
        %     z = zeros(length(this.z_opt), this.num_continuation_steps + 1);
        %     t = linspace(0, 1, this.num_continuation_steps + 1);

        %     u(:, 1) = this.u_opt;
        %     z(:, 1) = this.z_opt;

        %     for k = 1:this.num_continuation_steps
        %         Btheta_n = this.Apply_B(u(:, k), z(:, k), t(k));
        %         z_pert = -this.Apply_Parameterized_RS_Hessian_Inverse(Btheta_n, u(:, k), z(:, k), t(k));
        %         z_pred = z(:, k) + this.step_size * z_pert;
        %         u_pred = this.opt_prob_interface.State_Solve(z_pred);

        %         Jz_val = this.Gradient_J_z(u_pred, z_pred, t(k + 1));
        %         z(:, k + 1) = z_pred - this.Apply_Parameterized_RS_Hessian_Inverse(Jz_val, u_pred, z_pred, t(k + 1));
        %         u(:, k + 1) = this.opt_prob_interface.State_Solve(z(:, k + 1));
        %     end
        % end

        % function [u, z, beta] = Posterior_Update_Mean_beta(this)
        %     u = zeros(length(this.u_opt), this.num_continuation_steps + 1);
        %     z = zeros(length(this.z_opt), this.num_continuation_steps + 1);
        %     t = linspace(0, 1, this.num_continuation_steps + 1);
        %     beta = zeros(length(this.md_hessian_analysis.evals), this.num_continuation_steps + 1);

        %     u(:, 1) = this.u_opt;
        %     z(:, 1) = this.z_opt;

        %     for k = 1:this.num_continuation_steps
        %         Btheta_n = this.Apply_B_beta(u(:, k), beta(:, k), t(k));
        %         beta_pert = -this.Apply_Parameterized_RS_Hessian_Inverse_beta(Btheta_n, u(:, k), beta(:, k), t(k));
        %         beta(:, k + 1) = beta(:, k) + this.step_size * beta_pert;
        %         z(:, k + 1) = this.z_opt + this.md_hessian_analysis.evecs * beta(:, k + 1);
        %         u(:, k + 1) = this.opt_prob_interface.State_Solve(z(:, k + 1));
        %     end
        % end

        % function [z_out] = Apply_Parameterized_RS_Hessian_Inverse(this, z_in, u_n, z_n, t_n)
        %     z_out = 0 * z_in;
        %     for k = 1:size(z_in, 2)
        %         tol = 1.e-7;
        %         max_iter = length(z_n) + 5;
        %         [z_out(:, k), flag, relres, iter, resvec] = pcg(@(x)this.Apply_Parameterized_RS_Hessian(x, u_n, z_n, t_n), z_in(:, k), tol, max_iter);
        %         if flag ~= 0
        %             disp('CG did not converge');
        %         end
        %     end
        % end

        % function [beta_out] = Apply_Parameterized_RS_Hessian_Inverse_beta_noCG(this, beta_in, u_n, beta_n, t_n)
        %     I = eye(length(beta_n));
        %     H = I;
        %     for i = 1:size(beta_n, 1)
        %         H(:, i) = this.Apply_Parameterized_RS_Hessian_beta(I(:, i), u_n, beta_n, t_n);
        %     end
        %     beta_out = H \ beta_in;
        % end

    end

end
