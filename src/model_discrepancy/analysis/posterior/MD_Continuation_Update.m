classdef MD_Continuation_Update < handle

    properties
        md_post_sampling
        md_hessian_analysis
        opt_prob_interface
        u_opt
        z_opt
        num_continuation_steps
        r
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
            [u, z, beta] = this.Posterior_Update_Core(0);
        end

        function [u_ks, z_ks, beta_ks] = Posterior_Update_Samples(this)
            num_samples = this.md_post_sampling.post_data.num_samples;
            u_ks = zeros(length(this.u_opt), num_samples);
            z_ks = zeros(length(this.z_opt), num_samples);
            beta_ks = zeros(this.r, num_samples);

            for sample_idx = 1:num_samples
                [u_k, z_k, beta_k] = this.Posterior_Update_Core(sample_idx);
                u_ks(:, sample_idx) = u_k;
                z_ks(:, sample_idx) = z_k;
                beta_ks(:, sample_idx) = beta_k;
            end
        end

        % ------------------------------------------------------------
        % Core continuation driver (shared by mean and samples)
        % ------------------------------------------------------------

        function [u_k, z_k, beta_k] = Posterior_Update_Core(this, sample_idx)

            sen_op = MD_Continuation_Sensitivity_Operators(this.md_post_sampling, this.md_hessian_analysis);
            qn_prec = MD_Quasi_Newton_Preconditioner(this.md_hessian_analysis);
            beta_nom = zeros(this.r, 1);
            pt_cont = Pseudo_Time_Continuation(beta_nom, sen_op, qn_prec);
            theta_traj = MD_Discrepancy_Parameter_Trajectory(this.num_continuation_steps, sample_idx);
            beta_k_tmp = pt_cont.Pseudo_Time_Continuation_Forward_Euler(theta_traj);

            beta_k = beta_k_tmp(:, end);
            u_k = pt_cont.sen_op.current_u;
            z_k = pt_cont.sen_op.current_z;

        end

        % function [u_k, z_k, beta_k] = Posterior_Update_Core(this, sample_idx)
        %
        %     t = linspace(0, 1, this.num_continuation_steps + 1);
        %
        %     u_k = this.u_opt;
        %     beta_k = zeros(this.r, 1);
        %     z_k = this.z_opt;
        %
        %     for k = 1:this.num_continuation_steps
        %         % Predictive step for beta
        %         Btheta_n = this.Apply_B_beta(u_k, beta_k, t(k), sample_idx);
        %         beta_pert = -this.Apply_Parameterized_RS_Hessian_Inverse_beta(Btheta_n, u_k, beta_k, t(k), sample_idx);
        %         beta_pred = beta_k + (1 / this.num_continuation_steps) * beta_pert;
        %
        %         % Predictive step for z and u
        %         z_pred = this.z_opt + this.md_hessian_analysis.Apply_V(beta_pred);
        %         u_pred = this.opt_prob_interface.State_Solve(z_pred);
        %
        %         % Corrective step for beta
        %         grad = this.md_hessian_analysis.Apply_V_Transpose(this.Parameterized_RS_Objective_Gradient(u_pred, z_pred, t(k + 1), sample_idx));
        %         beta_k = beta_pred - this.Apply_Parameterized_RS_Hessian_Inverse_beta(grad, u_pred, beta_pred, t(k + 1), sample_idx);
        %
        %         % Update state
        %         z_k = this.z_opt + this.md_hessian_analysis.Apply_V(beta_k);
        %         u_k = this.opt_prob_interface.State_Solve(z_k);
        %     end
        % end

        % % ------------------------------------------------------------
        % % Continuation building blocks (shared)
        % % ------------------------------------------------------------
        %
        % function [Btheta_n] = Apply_B_beta(this, u_n, beta_n, t_n, sample_idx)
        %     z_n = this.z_opt + this.md_hessian_analysis.Apply_V(beta_n);
        %     Btheta_n = this.md_hessian_analysis.Apply_V_Transpose(this.Apply_B(u_n, z_n, t_n, sample_idx));
        % end
        %
        % function [beta_out] = Apply_Parameterized_RS_Hessian_Inverse_beta(this, beta_in, u_n, beta_n, t_n, sample_idx)
        %     beta_out = 0 * beta_in;
        %     for k = 1:size(beta_in, 2)
        %         tol = 1.e-7;
        %         max_iter = 2 * length(beta_n) + 10;
        %         [beta_out(:, k), flag, relres] = pcg( ...
        %                                              @(x)this.Apply_Parameterized_RS_Hessian_beta(x, u_n, beta_n, t_n, sample_idx), ...
        %                                              beta_in(:, k), tol, max_iter);
        %
        %         if flag ~= 0
        %             disp(['CG did not converge; flag: ', num2str(flag), ', relres: ', num2str(relres)]);
        %         end
        %     end
        % end
        %
        % function [beta_out] = Apply_Parameterized_RS_Hessian_beta(this, beta_in, u_n, beta_n, t_n, sample_idx)
        %     z_n = this.z_opt + this.md_hessian_analysis.Apply_V(beta_n);
        %     beta_out = this.md_hessian_analysis.Apply_V_Transpose(this.Apply_Parameterized_RS_Hessian(this.md_hessian_analysis.Apply_V(beta_in), u_n, z_n, t_n, sample_idx));
        % end
        %
        % function [grad] = Parameterized_RS_Objective_Gradient(this, u_n, z_n, t_n, sample_idx)
        %     disc_ops = this.Get_Discrepancy_Ops(sample_idx);
        %     delta = disc_ops.Eval(z_n, t_n);
        %     [~, grad_u, grad_z] = this.opt_prob_interface.Objective_Function(u_n + delta, z_n);
        %     z_tmp1 = disc_ops.Apply_z_Jacobian_Transpose(grad_u, z_n, t_n);
        %     z_tmp2 = this.opt_prob_interface.Apply_Solution_Operator_z_Jacobian_Transpose(grad_u, z_n);
        %     grad = grad_z + z_tmp1 + z_tmp2;
        % end
        %
        % function [Btheta_n] = Apply_B(this, u_n, z_n, t_n, sample_idx)
        %     disc_ops = this.Get_Discrepancy_Ops(sample_idx);
        %     delta = disc_ops.Eval(z_n, t_n);
        %
        %     u_tmp1 = disc_ops.Apply_theta_Jacobian(z_n);
        %     u_tmp2 = this.opt_prob_interface.Apply_Misfit_Hessian(u_tmp1, u_n + delta, z_n);
        %     z_tmp1 = this.opt_prob_interface.Apply_Solution_Operator_z_Jacobian_Transpose(u_tmp2, z_n);
        %
        %     z_tmp2 = disc_ops.Apply_z_Jacobian_Transpose(u_tmp2, z_n, t_n);
        %
        %     state_grad = this.opt_prob_interface.Misfit_Gradient(u_n + delta, z_n);
        %     z_tmp3 = disc_ops.Apply_z_theta_Hessian(state_grad, z_n);
        %
        %     Btheta_n = z_tmp1 + z_tmp2 + z_tmp3;
        % end
        %
        % function [z_out] = Apply_Parameterized_RS_Hessian(this, z_in, u_n, z_n, t_n, sample_idx)
        %     disc_ops = this.Get_Discrepancy_Ops(sample_idx);
        %     delta = disc_ops.Eval(z_n, t_n);
        %
        %     z_out = this.opt_prob_interface.Apply_RS_Hessian(z_in, z_n);
        %     z_out = z_out(:);
        %
        %     u_tmp1 = disc_ops.Apply_z_Jacobian(z_in, z_n, t_n);
        %     u_tmp2 = this.opt_prob_interface.Apply_Misfit_Hessian(u_tmp1, u_n + delta, z_n);
        %     z_out = z_out + this.opt_prob_interface.Apply_Solution_Operator_z_Jacobian_Transpose(u_tmp2, z_n);
        %
        %     z_out = z_out + disc_ops.Apply_z_Jacobian_Transpose(u_tmp2, z_n, t_n);
        %
        %     u_tmp3 = this.opt_prob_interface.Apply_Solution_Operator_z_Jacobian(z_in, z_n);
        %     u_tmp4 = this.opt_prob_interface.Apply_Misfit_Hessian(u_tmp3, u_n + delta, z_n);
        %     z_out = z_out + disc_ops.Apply_z_Jacobian_Transpose(u_tmp4, z_n, t_n);
        % end
        %
        % % ------------------------------------------------------------
        % % Discrepancy ops factory
        % % ------------------------------------------------------------
        %
        % function disc_ops = Get_Discrepancy_Ops(this, sample_idx)
        %     num_samples = this.md_post_sampling.post_data.num_samples;
        %     assert(sample_idx >= 0 && sample_idx <= num_samples && floor(sample_idx) == sample_idx, ...
        %            'sample_idx must be an integer in [0, num_samples].');
        %
        %     if sample_idx == 0 % Mean
        %         disc_ops.Eval          = @(z, t) t * this.md_post_sampling.Discrepancy_Evaluation_Mean(z);
        %         disc_ops.Apply_z_Jacobian       = @(z_in, z, t) t * this.md_post_sampling.Apply_Discrepancy_z_Jacobian_Mean(z_in);
        %         disc_ops.Apply_z_Jacobian_Transpose      = @(u_in, z, t) t * this.md_post_sampling.Apply_Discrepancy_z_Jacobian_Transpose_Mean(u_in);
        %         disc_ops.Apply_theta_Jacobian   = @(z) this.md_post_sampling.Discrepancy_Evaluation_Mean(z);
        %         disc_ops.Apply_z_theta_Hessian = @(u_in, z) this.md_post_sampling.Apply_Discrepancy_z_Jacobian_Transpose_Mean(u_in);
        %     else
        %         disc_ops.Eval          = @(z, t) t * this.md_post_sampling.Discrepancy_Evaluation_Sample(z, sample_idx);
        %         disc_ops.Apply_z_Jacobian       = @(z_in, z, t) t * this.md_post_sampling.Apply_Discrepancy_z_Jacobian_Sample(z_in, z, sample_idx);
        %         disc_ops.Apply_z_Jacobian_Transpose      = @(u_in, z, t) t * this.md_post_sampling.Apply_Discrepancy_z_Jacobian_Transpose_Sample(u_in, z, sample_idx);
        %         disc_ops.Apply_theta_Jacobian   = @(z) this.md_post_sampling.Discrepancy_Evaluation_Sample(z, sample_idx);
        %         disc_ops.Apply_z_theta_Hessian = @(u_in, z) this.md_post_sampling.Apply_Discrepancy_z_Jacobian_Transpose_Sample(u_in, z, sample_idx);
        %     end
        % end

        % ------------------------------------------------------------
        % Objective helpers (optional; kept for completeness)
        % ------------------------------------------------------------

        function [val, grad] = Parameterized_RS_Objective_beta(this, beta, sample_idx)
            z = this.z_opt + this.md_hessian_analysis.Apply_V(beta);
            [val, grad_z] = this.Parameterized_RS_Objective(z, sample_idx);
            grad = this.md_hessian_analysis.Apply_V_Transpose(grad_z);
        end

        function [val, grad] = Parameterized_RS_Objective(this, z, sample_idx)
            sen_op = MD_Continuation_Sensitivity_Operators(this.md_post_sampling, this.md_hessian_analysis);
            disc_ops = sen_op.Get_Discrepancy_Ops(sample_idx);
            u = this.opt_prob_interface.State_Solve(z);
            delta = disc_ops.Eval(z, 1);

            [val, grad_u, grad_z] = this.opt_prob_interface.Objective_Function(u + delta, z);
            z_tmp1 = disc_ops.Apply_z_Jacobian_Transpose(grad_u, z, 1);
            z_tmp2 = this.opt_prob_interface.Apply_Solution_Operator_z_Jacobian_Transpose(grad_u, z);
            grad = grad_z + z_tmp1 + z_tmp2;
        end

    end

end
