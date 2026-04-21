%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
