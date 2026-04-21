%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef MD_Update < handle

    properties
        post_sampling
        hessian_analysis
        opt_prob_interface
        u_opt
        z_opt
        state_grad
        state_grad_W_u_inv_state_grad
    end

    methods

        function this = MD_Update(post_sampling, hessian_analysis)
            arguments
                post_sampling MD_Posterior_Sampling
                hessian_analysis MD_Hessian_Analysis
            end
            this.post_sampling = post_sampling;
            this.hessian_analysis = hessian_analysis;
            this.opt_prob_interface = hessian_analysis.opt_prob_interface;
            this.u_opt = post_sampling.data_interface.u_opt;
            this.z_opt = post_sampling.data_interface.z_opt;

            this.state_grad = this.opt_prob_interface.Misfit_Gradient(this.u_opt, this.z_opt);
            this.state_grad_W_u_inv_state_grad = post_sampling.u_prior_interface.Apply_W_u_Inverse(this.state_grad)' * this.state_grad;
        end

        function [z_update_mean, z_update_samples] = Posterior_Update_Samples(this)
            z_update_mean = this.Posterior_Update_Mean();

            m = size(this.post_sampling.post_data.ui_hat{1}, 1);
            n = length(this.z_opt);

            Btheta_hat = zeros(n, this.post_sampling.post_data.num_samples);
            u_tmp1 = zeros(m, this.post_sampling.post_data.num_samples);
            for i = 1:this.post_sampling.post_data.N
                sgi = sum(this.post_sampling.post_data.g_vecs(:, i));
                coeff = sgi / sqrt(this.post_sampling.post_data.Mu(i, i));
                u_tmp1 = u_tmp1 + coeff * this.post_sampling.post_data.ui_hat{i};

                coeff = this.state_grad' * this.post_sampling.post_data.ui_hat{i};
                M_z_W_z_inv_M_z_yi = this.post_sampling.post_data.Mz_Wz_inv_Mz_Z * this.post_sampling.post_data.g_vecs(:, i) - sgi * this.post_sampling.post_data.Mz_Wz_inv_Mz_z_opt;
                Btheta_hat = Btheta_hat + (1 / sqrt(this.post_sampling.post_data.Mu(i, i))) * M_z_W_z_inv_M_z_yi * coeff;
            end
            tmp = this.opt_prob_interface.Apply_Misfit_Hessian(u_tmp1, this.u_opt, this.z_opt);
            Btheta_hat = Btheta_hat + this.opt_prob_interface.Apply_Solution_Operator_z_Jacobian_Transpose(tmp, this.z_opt);
            Btheta_hat = sqrt(this.post_sampling.post_data.alpha_d) * Btheta_hat;

            tmp1 = this.post_sampling.post_data.Zc' * this.post_sampling.post_data.Mz_z_breve;
            tmp2 = linsolve(this.post_sampling.post_data.Zc_Mz_Wz_inv_Mz_Zc, tmp1);
            tmp3 = this.post_sampling.post_data.Mz_Wz_inv_Mz_Zc * tmp2;
            coeff = sqrt(this.state_grad_W_u_inv_state_grad);
            Btheta_breve = coeff * (this.post_sampling.post_data.Mz_z_breve - tmp3);

            z_update_samples = z_update_mean - this.hessian_analysis.Apply_RS_Hessian_Inverse(Btheta_hat + Btheta_breve, this.z_opt);
        end

        function [z_update_mean] = Posterior_Update_Mean(this)

            N = this.post_sampling.post_data.N;
            u = 0 * this.u_opt;
            for ell = 1:N
                u = u + this.post_sampling.post_data.u_ell(:, ell);
                for i = 1:N
                    u = u - this.post_sampling.post_data.b_i_ell(i, ell) * sum(this.post_sampling.post_data.g_vecs(:, i)) * this.post_sampling.post_data.u_i_ell{i}(:, ell);
                end
            end
            u = u  + this.post_sampling.post_data.alpha_d * this.post_sampling.data_interface.data_shift;
            tmp1 = this.opt_prob_interface.Apply_Misfit_Hessian(u, this.u_opt, this.z_opt);
            z_tmp = this.opt_prob_interface.Apply_Solution_Operator_z_Jacobian_Transpose(tmp1, this.z_opt);

            for ell = 1:N
                z_tmp = z_tmp + (this.state_grad' * this.post_sampling.post_data.u_ell(:, ell)) * (this.post_sampling.post_data.Mz_Wz_inv_Mz_Z(:, ell) - this.post_sampling.post_data.Mz_Wz_inv_Mz_z_opt);
                for i = 1:N
                    coeff = this.post_sampling.post_data.b_i_ell(i, ell) * (this.state_grad' * this.post_sampling.post_data.u_i_ell{i}(:, ell));
                    vec = this.post_sampling.post_data.Mz_Wz_inv_Mz_Z * this.post_sampling.post_data.g_vecs(:, i) - sum(this.post_sampling.post_data.g_vecs(:, i)) * this.post_sampling.post_data.Mz_Wz_inv_Mz_z_opt;
                    z_tmp = z_tmp - coeff * vec;
                end
            end

            z_tmp = (1 / this.post_sampling.post_data.alpha_d) * z_tmp;
            z_pert = -this.hessian_analysis.Apply_RS_Hessian_Inverse(z_tmp, this.z_opt);
            z_update_mean = this.z_opt + z_pert;
        end

    end

end
