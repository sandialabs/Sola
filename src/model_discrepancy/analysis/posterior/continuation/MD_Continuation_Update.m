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
    end

    methods

        function this = MD_Continuation_Update(post_sampling, hessian_analysis, num_continuation_steps)
            arguments
                post_sampling MD_Posterior_Sampling
                hessian_analysis MD_Hessian_Analysis
                num_continuation_steps (1, 1) {mustBeNumeric}
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
        end

        % ------------------------------------------------------------
        % High-level functions
        % ------------------------------------------------------------

        function [u, z, beta] = Posterior_Update_Mean(this)
            [u, z, beta] = this.Posterior_Update_Core(0);
        end

        function [u_ks, z_ks, beta_ks] = Posterior_Update_Samples(this)
            num_samples = this.post_sampling.post_data.num_samples;
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

            % NOTE: All continuation classes prefixed with `MD_` update the latent variable beta, 
            % although denoted by parent classes as z. This is not to be confused with the solution
            % space optimization variable z in model discrepancy calibration.
            sen_op = MD_Continuation_Sensitivity_Operators(this.post_sampling, this.hessian_analysis);
            qn_prec = MD_Quasi_Newton_Preconditioner(this.hessian_analysis);
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
            beta = beta(:);
            sen_op = MD_Continuation_Sensitivity_Operators(this.post_sampling, this.hessian_analysis);
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