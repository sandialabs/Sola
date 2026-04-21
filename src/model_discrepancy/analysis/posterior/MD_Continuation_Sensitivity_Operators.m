%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef MD_Continuation_Sensitivity_Operators < Sensitivity_Operators

    properties
        md_post_sampling
        post_data
        md_hessian_analysis
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
    end

    methods (Access = public)

        function [grad, val] = Gradient(this, z, theta_traj, time_index)

            this.State_Evaluation(z, theta_traj, time_index);

            delta = this.current_disc_ops.Eval(this.current_z, this.current_t);
            [val, grad_u, grad_z] = this.opt_prob_interface.Objective_Function(this.current_u + delta, this.current_z);
            z_tmp1 = this.current_disc_ops.Apply_z_Jacobian_Transpose(grad_u, this.current_z, this.current_t);
            z_tmp2 = this.opt_prob_interface.Apply_Solution_Operator_z_Jacobian_Transpose(grad_u, this.current_z);
            grad = grad_z + z_tmp1 + z_tmp2;

            grad = this.md_hessian_analysis.Apply_V_Transpose(grad);
        end

        function [z_out] = Apply_Hessian(this, z_in, z, theta_traj, time_index)

            this.State_Evaluation(z, theta_traj, time_index);
            V_z_in = this.md_hessian_analysis.Apply_V(z_in);
            delta = this.current_disc_ops.Eval(this.current_z, this.current_t);

            z_out = this.opt_prob_interface.Apply_RS_Hessian(V_z_in, this.current_z);

            u_tmp1 = this.current_disc_ops.Apply_z_Jacobian(V_z_in,  this.current_z,  this.current_t);
            u_tmp2 = this.opt_prob_interface.Apply_Misfit_Hessian(u_tmp1,  this.current_u + delta,  this.current_z);
            z_out = z_out + this.opt_prob_interface.Apply_Solution_Operator_z_Jacobian_Transpose(u_tmp2,  this.current_z);

            z_out = z_out + this.current_disc_ops.Apply_z_Jacobian_Transpose(u_tmp2,  this.current_z,  this.current_t);

            u_tmp3 = this.opt_prob_interface.Apply_Solution_Operator_z_Jacobian(V_z_in,  this.current_z);
            u_tmp4 = this.opt_prob_interface.Apply_Misfit_Hessian(u_tmp3,  this.current_u + delta,  this.current_z);
            z_out = z_out + this.current_disc_ops.Apply_z_Jacobian_Transpose(u_tmp4,  this.current_z,  this.current_t);

            z_out = this.md_hessian_analysis.Apply_V_Transpose(z_out);

        end

        function [z_out] = Apply_B(this, z, theta_traj, time_index)

            this.State_Evaluation(z, theta_traj, time_index);
            delta = this.current_disc_ops.Eval(this.current_z, this.current_t);

            u_tmp1 = this.current_disc_ops.Apply_theta_Jacobian(this.current_z);
            u_tmp2 = this.opt_prob_interface.Apply_Misfit_Hessian(u_tmp1, this.current_u + delta, this.current_z);
            z_tmp1 = this.opt_prob_interface.Apply_Solution_Operator_z_Jacobian_Transpose(u_tmp2, this.current_z);

            z_tmp2 = this.current_disc_ops.Apply_z_Jacobian_Transpose(u_tmp2, this.current_z, this.current_t);

            state_grad = this.opt_prob_interface.Misfit_Gradient(this.current_u + delta, this.current_z);
            z_tmp3 = this.current_disc_ops.Apply_z_theta_Hessian(state_grad, this.current_z);

            z_out = this.md_hessian_analysis.Apply_V_Transpose(z_tmp1 + z_tmp2 + z_tmp3);

        end

    end

    methods

        function this = MD_Continuation_Sensitivity_Operators(md_post_sampling, md_hessian_analysis)
            arguments
                md_post_sampling MD_Posterior_Sampling
                md_hessian_analysis MD_Hessian_Analysis
            end

            this.md_post_sampling = md_post_sampling;
            this.post_data = md_post_sampling.post_data;
            this.md_hessian_analysis = md_hessian_analysis;
            this.opt_prob_interface = md_hessian_analysis.opt_prob_interface;
            this.data_interface = md_post_sampling.data_interface;
            this.z_prior_interface = md_post_sampling.z_prior_interface;
            this.u_opt = md_post_sampling.data_interface.u_opt;
            this.z_opt = md_post_sampling.data_interface.z_opt;

            this.current_t = inf;
            this.current_beta = inf;

            this.Mz_Wz_inv_Mz_Z_minus_z_opt = this.post_data.Mz_Wz_inv_Mz_Z - this.post_data.Mz_Wz_inv_Mz_z_opt;
            this.Mz_Wz_inv_Mz_yi = 0 * this.Mz_Wz_inv_Mz_Z_minus_z_opt;
            this.si = zeros(1, this.post_data.N);
            for i = 1:this.post_data.N
                this.Mz_Wz_inv_Mz_yi(:, i) = this.post_data.Mz_Wz_inv_Mz_Z * this.post_data.g_vecs(:, i) - ...
                    sum(this.post_data.g_vecs(:, i)) * this.post_data.Mz_Wz_inv_Mz_z_opt;
                this.si(i) = sum(this.post_data.g_vecs(:, i)) - this.z_opt' * this.Mz_Wz_inv_Mz_yi(:, i);
            end
        end

        function [] = State_Evaluation(this, beta, theta_traj, time_index)
            t = theta_traj.Get_Time(time_index);
            if max(abs(t - this.current_t), norm(beta - this.current_beta)) > 1.e-15
                this.current_t = t;
                this.current_beta = beta;
                this.current_z = this.z_opt + this.md_hessian_analysis.Apply_V(beta);
                this.current_u = this.opt_prob_interface.State_Solve(this.current_z);
                this.current_disc_ops = this.Get_Discrepancy_Ops(theta_traj.Get_Sample_Index());
            end
        end

        function [disc_ops] = Get_Discrepancy_Ops(this, sample_idx)
            num_samples = this.post_data.num_samples;
            assert(sample_idx >= 0 && sample_idx <= num_samples && floor(sample_idx) == sample_idx, ...
                   'sample_idx must be an integer in [0, num_samples].');

            if sample_idx == 0 % Mean
                disc_ops.Eval          = @(z, t) t * this.Discrepancy_Evaluation_Mean(z);
                disc_ops.Apply_z_Jacobian       = @(z_in, z, t) t * this.Apply_Discrepancy_z_Jacobian_Mean(z_in);
                disc_ops.Apply_z_Jacobian_Transpose      = @(u_in, z, t) t * this.Apply_Discrepancy_z_Jacobian_Transpose_Mean(u_in);
                disc_ops.Apply_theta_Jacobian   = @(z) this.Discrepancy_Evaluation_Mean(z);
                disc_ops.Apply_z_theta_Hessian = @(u_in, z) this.Apply_Discrepancy_z_Jacobian_Transpose_Mean(u_in);
            else
                disc_ops.Eval          = @(z, t) t * this.Discrepancy_Evaluation_Sample(z, sample_idx);
                disc_ops.Apply_z_Jacobian       = @(z_in, z, t) t * this.Apply_Discrepancy_z_Jacobian_Sample(z_in, z, sample_idx);
                disc_ops.Apply_z_Jacobian_Transpose      = @(u_in, z, t) t * this.Apply_Discrepancy_z_Jacobian_Transpose_Sample(u_in, z, sample_idx);
                disc_ops.Apply_theta_Jacobian   = @(z) this.Discrepancy_Evaluation_Sample(z, sample_idx);
                disc_ops.Apply_z_theta_Hessian = @(u_in, z) this.Apply_Discrepancy_z_Jacobian_Transpose_Sample(u_in, z, sample_idx);
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
        % Discrepancy kernels: Sample
        % ------------------------------------------------------------

        function [u_out] = Discrepancy_Evaluation_Sample(this, z, sample_idx)
            u_out_mean = this.Discrepancy_Evaluation_Mean(z);
            dz = z - this.z_opt;

            Mz_dz = this.z_prior_interface.Apply_M_z(dz);
            Wz_inv_Mz_dz = this.z_prior_interface.Apply_W_z_Inverse(Mz_dz);

            delta_sample = zeros(size(u_out_mean));
            for i = 1:this.post_data.N
                sgi = sum(this.post_data.g_vecs(:, i));
                coeff = (1 / sqrt(this.post_data.Mu(i, i))) * (sgi + this.Mz_Wz_inv_Mz_yi(:, i)' * dz);
                delta_sample = delta_sample + coeff * this.post_data.ui_hat{i}(:, sample_idx);
            end
            delta_sample = sqrt(this.post_data.alpha_d) * delta_sample;

            tmp = Mz_dz' * Wz_inv_Mz_dz - ...
                Wz_inv_Mz_dz' * this.post_data.Mz_Zc * linsolve( ...
                                                                this.post_data.Zc_Mz_Wz_inv_Mz_Zc, ...
                                                                this.post_data.Mz_Zc' * Wz_inv_Mz_dz);

            if tmp < -1.e-11
                disp('Error in Posterior Discrepancy Sample: delta breve coeff < 0');
            end

            breve_coeff = sqrt(abs(tmp));
            delta_sample = delta_sample + breve_coeff * this.post_data.u_breve(:, sample_idx);

            u_out = u_out_mean + delta_sample;
        end

        function [u_out] = Apply_Discrepancy_z_Jacobian_Sample(this, z_in, z, sample_idx)
            % Note: z is needed since sampling is nonlinear in z via gamma(z) in delta_breve
            u_out_mean = this.Apply_Discrepancy_z_Jacobian_Mean(z_in);

            u = zeros(size(u_out_mean));
            for i = 1:this.post_data.N
                coeff = (1 / sqrt(this.post_data.Mu(i, i))) * (this.Mz_Wz_inv_Mz_yi(:, i)' * z_in);
                u = u + coeff * this.post_data.ui_hat{i}(:, sample_idx);
            end
            u = sqrt(this.post_data.alpha_d) * u;

            Mz_dz = this.z_prior_interface.Apply_M_z(z - this.z_opt);
            Wz_inv_Mz_dz = this.z_prior_interface.Apply_W_z_Inverse(Mz_dz);

            tmp_rhs = Wz_inv_Mz_dz - this.post_data.Wz_inv_Mz_Zc * linsolve( ...
                                                                            this.post_data.Zc_Mz_Wz_inv_Mz_Zc, ...
                                                                            this.post_data.Mz_Zc' * Wz_inv_Mz_dz);

            tmp = Mz_dz' * tmp_rhs;
            if tmp < -1.e-11
                disp('Error in Posterior Discrepancy Samples: delta breve coeff < 0');
            end

            Mz_z_in = this.z_prior_interface.Apply_M_z(z_in);
            denom = sqrt(abs(tmp) + (1e-15)^2); % tiny regularization
            breve_coeff_deriv = (Mz_z_in' * tmp_rhs) / denom;
            u = u + breve_coeff_deriv * this.post_data.u_breve(:, sample_idx);

            u_out = u_out_mean + u;
        end

        function [z_out] = Apply_Discrepancy_z_Jacobian_Transpose_Sample(this, u_in, z, sample_idx)
            % Note: z is needed since sampling is nonlinear in z via gamma(z) in delta_breve
            z_out_mean = this.Apply_Discrepancy_z_Jacobian_Transpose_Mean(u_in);

            z_out_sample = zeros(size(z_out_mean));
            for i = 1:this.post_data.N
                ui_hat_idx = this.post_data.ui_hat{i}(:, sample_idx);
                coeff = (1 / sqrt(this.post_data.Mu(i, i))) * (ui_hat_idx' * u_in);
                z_out_sample = z_out_sample + coeff * this.Mz_Wz_inv_Mz_yi(:, i);
            end
            z_out_sample = sqrt(this.post_data.alpha_d) * z_out_sample;

            Mz_dz = this.z_prior_interface.Apply_M_z(z - this.z_opt);
            Wz_inv_Mz_dz = this.z_prior_interface.Apply_W_z_Inverse(Mz_dz);

            tmp_rhs = Wz_inv_Mz_dz - this.post_data.Wz_inv_Mz_Zc * linsolve( ...
                                                                            this.post_data.Zc_Mz_Wz_inv_Mz_Zc, ...
                                                                            this.post_data.Mz_Zc' * Wz_inv_Mz_dz);

            tmp = Mz_dz' * tmp_rhs;
            if tmp < -1.e-11
                disp('Error in Posterior Discrepancy Samples: delta breve coeff < 0');
            end

            denom = sqrt(abs(tmp) + (1e-15)^2); % tiny regularization
            breve_coeff_grad = this.z_prior_interface.Apply_M_z(tmp_rhs) / denom;

            u_breve_idx = this.post_data.u_breve(:, sample_idx);
            z_out_sample = z_out_sample + (u_breve_idx' * u_in) * breve_coeff_grad;

            z_out = z_out_mean + z_out_sample;
        end

    end
end
