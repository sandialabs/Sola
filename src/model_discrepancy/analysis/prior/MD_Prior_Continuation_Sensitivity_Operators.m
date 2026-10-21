%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef MD_Prior_Continuation_Sensitivity_Operators < Sensitivity_Operators

    % Sensitivity operators for continuation of a prior discrepancy sample.
    %
    % For sample k, the discrepancy at pseudo-time t is
    %
    %   delta_k(beta, t) = t * (data_shift + delta0_k + X_k * beta),
    %
    % where delta0_k ~ N(0, W_u^{-1}) and X_k ~ MN(0, W_u^{-1}, Sigma_beta)
    % are fixed realizations (see MD_Prior_Continuation_Sampling).  At t = 0
    % the discrepancy vanishes and beta = 0 is optimal.

    properties
        hessian_analysis
        opt_prob_interface
        data_shift
        delta0
        X

        current_t
        current_u
        current_beta
        current_z
        current_sample_idx
        current_disc_ops
        z_opt
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

            % Corrections to the low-fidelity reduced Hessian (see
            % MD_Continuation_Sensitivity_Operators::Apply_Hessian).
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

        % delta0: m x num_samples, X: cell(num_samples, 1) of m x r matrices
        function this = MD_Prior_Continuation_Sensitivity_Operators(hessian_analysis, z_opt, data_shift, delta0, X)
            arguments
                hessian_analysis MD_Hessian_Analysis
                z_opt (:, 1) {mustBeNumeric}
                data_shift (:, 1) {mustBeNumeric}
                delta0 {mustBeNumeric}
                X cell
            end
            this.hessian_analysis = hessian_analysis;
            this.opt_prob_interface = hessian_analysis.opt_prob_interface;
            this.data_shift = data_shift;
            this.delta0 = delta0;
            this.X = X;
            this.z_opt = z_opt;

            this.current_t = inf;
            this.current_beta = inf;
            this.current_sample_idx = inf;
        end

        function [] = State_Evaluation(this, beta, theta_traj, time_index)
            beta = beta(:);

            t = theta_traj.Get_Time(time_index);
            sample_idx = theta_traj.Get_Sample_Index();

            if max([abs(t - this.current_t), norm(beta - this.current_beta), abs(sample_idx - this.current_sample_idx)]) > 1.e-15
                this.current_t = t;
                this.current_beta = beta;
                this.current_sample_idx = sample_idx;
                this.current_z = this.z_opt + this.hessian_analysis.Apply_V(beta);
                this.current_u = this.opt_prob_interface.State_Solve(this.current_z);
                this.current_disc_ops = this.Get_Discrepancy_Ops(sample_idx);
            end
        end

        function [disc_ops] = Get_Discrepancy_Ops(this, sample_idx)
            assert(sample_idx >= 1 && sample_idx <= size(this.delta0, 2) && floor(sample_idx) == sample_idx, ...
                'sample_idx must be an integer in [1, num_samples].');

            X_k = this.X{sample_idx};
            d0 = this.data_shift + this.delta0(:, sample_idx);
            disc_ops.Eval = @(beta, t) t * (d0 + X_k * beta);
            disc_ops.Apply_Beta_Jacobian = @(beta_in, t) t * (X_k * beta_in);
            disc_ops.Apply_Beta_Jacobian_Transpose = @(u_in, t) t * (X_k' * u_in);
            disc_ops.Apply_Theta_Jacobian = @(beta) d0 + X_k * beta;
            disc_ops.Apply_Beta_Theta_Hessian = @(u_in) X_k' * u_in;
        end

    end

end
