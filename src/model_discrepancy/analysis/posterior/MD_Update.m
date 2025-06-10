classdef MD_Update < handle

    properties
        md_post_sampling
        md_hessian_analysis
        opt_prob_interface
        u_init
        z_init
        state_grad
        state_grad_W_u_inv_state_grad
    end

    methods

        function this = MD_Update(md_post_sampling, md_hessian_analysis)
            arguments
                md_post_sampling MD_Posterior_Sampling
                md_hessian_analysis MD_Hessian_Analysis
            end
            this.md_post_sampling = md_post_sampling;
            this.md_hessian_analysis = md_hessian_analysis;
            this.opt_prob_interface = md_hessian_analysis.opt_prob_interface;
            this.u_init = md_post_sampling.data_interface.u_init;
            this.z_init = md_post_sampling.data_interface.z_init;

            this.state_grad = this.opt_prob_interface.Misfit_Gradient(this.u_init, this.z_init);
            this.state_grad_W_u_inv_state_grad = md_post_sampling.u_prior_interface.Apply_W_u_Inverse(this.state_grad)' * this.state_grad;
        end

        function [z_update_mean, z_update_samples] = Posterior_Update_Samples(this)
            z_update_mean = this.Posterior_Update_Mean();

            m = size(this.md_post_sampling.post_data.ui_hat{1}, 1);
            n = length(this.z_init);

            Btheta_hat = zeros(n, this.md_post_sampling.post_data.num_samples);
            u_tmp1 = zeros(m, this.md_post_sampling.post_data.num_samples);
            for i = 1:this.md_post_sampling.post_data.N
                sgi = sum(this.md_post_sampling.post_data.g_vecs(:, i));
                coeff = sgi / sqrt(this.md_post_sampling.post_data.Mu(i, i));
                u_tmp1 = u_tmp1 + coeff * this.md_post_sampling.post_data.ui_hat{i};

                coeff = this.state_grad' * this.md_post_sampling.post_data.ui_hat{i};
                M_z_W_z_inv_M_z_yi = this.md_post_sampling.post_data.Mz_Wz_inv_Mz_Z * this.md_post_sampling.post_data.g_vecs(:, i) - sgi * this.md_post_sampling.post_data.Mz_Wz_inv_Mz_z_opt;
                Btheta_hat = Btheta_hat + (1 / sqrt(this.md_post_sampling.post_data.Mu(i, i))) * M_z_W_z_inv_M_z_yi * coeff;
            end
            tmp = this.opt_prob_interface.Apply_Misfit_Hessian(u_tmp1, this.u_init, this.z_init);
            Btheta_hat = Btheta_hat + this.opt_prob_interface.Apply_Solution_Operator_z_Jacobian_Transpose(tmp, this.z_init);
            Btheta_hat = sqrt(this.md_post_sampling.post_data.alpha_d) * Btheta_hat;

            tmp1 = this.md_post_sampling.post_data.Zc' * this.md_post_sampling.post_data.Mz_z_breve;
            tmp2 = linsolve(this.md_post_sampling.post_data.Zc_Mz_Wz_inv_Mz_Zc, tmp1);
            tmp3 = this.md_post_sampling.post_data.Mz_Wz_inv_Mz_Zc * tmp2;
            coeff = sqrt(this.state_grad_W_u_inv_state_grad);
            Btheta_breve = coeff * (this.md_post_sampling.post_data.Mz_z_breve - tmp3);

            z_update_samples = z_update_mean - this.md_hessian_analysis.Apply_RS_Hessian_Inverse(Btheta_hat + Btheta_breve, this.z_init);
        end

        function [z_update_mean] = Posterior_Update_Mean(this)

            N = this.md_post_sampling.post_data.N;
            diff_z_opt = this.z_init - this.md_post_sampling.data_interface.z_opt; % set to zero to reproduce initial bug
            % diff_z_opt = 0*diff_z_opt; % to reproduce initial bug
            u = 0 * this.u_init;
            for ell = 1:N
                coeff = 1 + diff_z_opt' * (this.md_post_sampling.post_data.Mz_Wz_inv_Mz_Z(:, ell) - this.md_post_sampling.post_data.Mz_Wz_inv_Mz_z_opt);
                u = u + coeff * this.md_post_sampling.post_data.u_ell(:, ell);
                for i = 1:N
                    vec = this.md_post_sampling.post_data.Mz_Wz_inv_Mz_Z * this.md_post_sampling.post_data.g_vecs(:, i) - sum(this.md_post_sampling.post_data.g_vecs(:, i)) * this.md_post_sampling.post_data.Mz_Wz_inv_Mz_z_opt;
                    coeff = sum(this.md_post_sampling.post_data.g_vecs(:, i)) + diff_z_opt' * vec;
                    u = u - coeff * this.md_post_sampling.post_data.b_i_ell(i, ell) * this.md_post_sampling.post_data.u_i_ell{i}(:, ell);
                end
            end

            u = u + this.md_post_sampling.post_data.alpha_d * this.md_post_sampling.data_interface.data_shift;
            tmp1 = this.opt_prob_interface.Apply_Misfit_Hessian(u, this.u_init, this.z_init);
            z_tmp = this.opt_prob_interface.Apply_Solution_Operator_z_Jacobian_Transpose(tmp1, this.z_init);

            for ell = 1:N
                z_tmp = z_tmp + (this.state_grad' * this.md_post_sampling.post_data.u_ell(:, ell)) * (this.md_post_sampling.post_data.Mz_Wz_inv_Mz_Z(:, ell) - this.md_post_sampling.post_data.Mz_Wz_inv_Mz_z_opt);
                for i = 1:N
                    coeff = this.md_post_sampling.post_data.b_i_ell(i, ell) * (this.state_grad' * this.md_post_sampling.post_data.u_i_ell{i}(:, ell));
                    vec = this.md_post_sampling.post_data.Mz_Wz_inv_Mz_Z * this.md_post_sampling.post_data.g_vecs(:, i) - sum(this.md_post_sampling.post_data.g_vecs(:, i)) * this.md_post_sampling.post_data.Mz_Wz_inv_Mz_z_opt;
                    z_tmp = z_tmp - coeff * vec;
                end
            end

            z_tmp = (1 / this.md_post_sampling.post_data.alpha_d) * z_tmp;
            z_pert = -this.md_hessian_analysis.Apply_RS_Hessian_Inverse(z_tmp, this.z_init);
            z_update_mean = this.z_init + z_pert;
        end

    end

end
