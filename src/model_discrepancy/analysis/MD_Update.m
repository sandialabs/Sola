classdef MD_Update < handle

    properties
        opt_prob_interface
        data_interface
        u_prior_interface
        z_prior_interface
        md_hessian_analysis
        post_data
        u_opt
        z_opt
    end

    methods

        function this = MD_Update(opt_prob_interface, data_interface, u_prior_interface, z_prior_interface, md_hessian_analysis)
            this.opt_prob_interface = opt_prob_interface;
            this.data_interface = data_interface;
            this.u_prior_interface = u_prior_interface;
            this.z_prior_interface = z_prior_interface;
            this.md_hessian_analysis = md_hessian_analysis;
            this.post_data = MD_Bayes_Posterior_Data();
            this.u_opt = data_interface.u_opt;
            this.z_opt = data_interface.z_opt;
        end

        function [] = Compute_Posterior_Data(this, alpha_d, num_samples)
            this.post_data.Compute_Posterior_Data(this.opt_prob_interface, this.data_interface, this.u_prior_interface, this.z_prior_interface, alpha_d, this.u_opt, this.z_opt, num_samples);
        end

        function [mean_error, sample_error, delta_mean, delta_samples] = Compute_Discrepancy_Fit_Error(this)

            [delta_mean, delta_samples] = this.Posterior_Discrepancy_Samples(this.post_data.Z);
            mean_error = zeros(this.post_data.N, 1);
            sample_error = zeros(this.post_data.N, this.post_data.num_samples);
            for i = 1:this.post_data.N
                diff = delta_mean{i} - this.post_data.D(:, i);
                normalization = sqrt(this.post_data.D(:, i)' * this.u_prior_interface.Apply_M_u(this.post_data.D(:, i)));
                mean_error(i) = sqrt(this.u_prior_interface.Apply_M_u(diff)' * diff) / normalization;
                for k = 1:this.post_data.num_samples
                    diff = delta_samples{i}(:, k) - this.post_data.D(:, i);
                    sample_error(i, k) = sqrt(this.u_prior_interface.Apply_M_u(diff)' * diff) / normalization;
                end
            end

            for i = 1:this.post_data.N
                figure;
                hold on;
                histogram(sample_error(i, :));
                plot([mean_error(i), mean_error(i)], [0, max(get(gca, 'YLim'))], 'LineWidth', 5);
                title(['Posterior discrepancy relative error at $z_\ell$, $\ell=$', num2str(i)], 'Interpreter', 'latex');
                legend({'Posterior Samples', 'Posterior Mean'}, 'Location', 'best');
                set(gca, 'fontsize', 18);
            end
        end

        function [] = Compute_Discrepancy_Extrapolation_Variabilty(this, Z_test)

            [delta_mean, delta_samples] = this.Posterior_Discrepancy_Samples(Z_test);
            N = size(Z_test, 2);
            sample_variability = zeros(N, this.post_data.num_samples);
            for i = 1:N
                normalization = sqrt(delta_mean{i}' * this.u_prior_interface.Apply_M_u(delta_mean{i}));
                for k = 1:this.post_data.num_samples
                    diff = delta_samples{i}(:, k) - delta_mean{i};
                    sample_variability(i, k) = sqrt(this.u_prior_interface.Apply_M_u(diff)' * diff) / normalization;
                end
            end

            for i = 1:N
                figure;
                hold on;
                histogram(sample_variability(i, :));
                title(['Posterior discrepancy sample variability at $z_{test,i}$, $i=$', num2str(i)], 'Interpreter', 'latex');
                set(gca, 'fontsize', 18);
            end
        end

        function [delta_mean, delta_samples] = Posterior_Discrepancy_Samples(this, z)
            m = size(this.post_data.u_ell, 1);
            p = size(z, 2);
            delta_mean = cell(p, 1);
            delta_samples = cell(p, 1);

            Zc = this.post_data.Z(:, 2:end) - this.z_opt;
            Zc_W_z_Inv_Zc = Zc' * (this.post_data.W_z_inv_Z(:, 2:end) - this.post_data.W_z_inv_z_opt);

            for k = 1:p
                dz = z(:, k) - this.z_opt;

                coeffs = 1 + (this.post_data.W_z_inv_Z - this.post_data.W_z_inv_z_opt)' * dz;
                delta_mean_k = this.post_data.u_ell * coeffs;

                delta_samples_k = zeros(m, this.post_data.num_samples);
                for i = 1:this.post_data.N
                    sgi = sum(this.post_data.g_vecs(:, i));
                    W_z_Inv_yi = this.post_data.W_z_inv_Z * this.post_data.g_vecs(:, i) - sgi * this.post_data.W_z_inv_z_opt;

                    coeffs = this.post_data.b_i_ell(i, :)' * (sgi + W_z_Inv_yi' * dz);
                    delta_mean_k = delta_mean_k - this.post_data.u_i_ell{i} * coeffs;

                    coeff = (1 / sqrt(this.post_data.Mu(i, i))) * (sgi + W_z_Inv_yi' * dz);
                    delta_samples_k = delta_samples_k + coeff * this.post_data.ui_hat{i};
                end
                delta_mean_k = (1 / this.post_data.alpha_d) * delta_mean_k;
                delta_samples_k = sqrt(this.post_data.alpha_d) * delta_samples_k;

                W_z_Inv_dz = this.z_prior_interface.Apply_W_z_Inverse(dz);
                tmp = dz' * W_z_Inv_dz - W_z_Inv_dz' * Zc * linsolve(Zc_W_z_Inv_Zc, Zc' * W_z_Inv_dz);
                if tmp < -1.e-13
                    disp('Error in Posterior Discrepancy Samples');
                end
                breve_coeff = sqrt(abs(tmp));
                delta_samples_k = delta_samples_k + breve_coeff * this.post_data.u_breve;

                delta_mean{k} = delta_mean_k;
                delta_samples{k} = delta_samples_k + delta_mean_k;

            end
        end

        function [z_update_mean, z_update_samples] = Posterior_Update_Samples(this)
            z_update_mean = this.Posterior_Update_Mean();

            m = size(this.post_data.ui_hat{1}, 1);
            n = length(this.z_opt);

            Btheta_hat = zeros(n, this.post_data.num_samples);
            u_tmp1 = zeros(m, this.post_data.num_samples);
            for i = 1:this.post_data.N
                sgi = sum(this.post_data.g_vecs(:, i));
                coeff = sgi / sqrt(this.post_data.Mu(i, i));
                u_tmp1 = u_tmp1 + coeff * this.post_data.ui_hat{i};

                coeff = this.post_data.state_grad' * this.post_data.ui_hat{i};
                W_z_Inv_yi = this.post_data.W_z_inv_Z * this.post_data.g_vecs(:, i) - sgi * this.post_data.W_z_inv_z_opt;
                Btheta_hat = Btheta_hat + (1 / sqrt(this.post_data.Mu(i, i))) * W_z_Inv_yi * coeff;
            end
            tmp = this.opt_prob_interface.Apply_Misfit_Hessian(u_tmp1, this.u_opt, this.z_opt);
            Btheta_hat = Btheta_hat + this.opt_prob_interface.Apply_Solution_Operator_z_Jacobian_Transpose(tmp, this.z_opt);
            Btheta_hat = sqrt(this.post_data.alpha_d) * Btheta_hat;

            Zc = this.post_data.Z(:, 2:end) - this.z_opt;
            W_z_Inv_Zc = this.post_data.W_z_inv_Z(:, 2:end) - this.post_data.W_z_inv_z_opt;
            Zc_W_z_Inv_Zc = Zc' * W_z_Inv_Zc;
            tmp1 = Zc' * this.post_data.zbreve;
            tmp2 = linsolve(Zc_W_z_Inv_Zc, tmp1);
            tmp3 = W_z_Inv_Zc * tmp2;
            coeff = sqrt(this.post_data.state_grad_W_u_inv_state_grad);
            Btheta_breve = coeff * (this.post_data.zbreve - tmp3);

            z_update_samples = z_update_mean - this.md_hessian_analysis.Apply_RS_Hessian_Inverse(Btheta_hat + Btheta_breve, this.z_opt);
        end

        function [z_update_mean] = Posterior_Update_Mean(this)

            N = this.post_data.N;
            u = 0 * this.u_opt;
            for ell = 1:N
                u = u + this.post_data.u_ell(:, ell);
                for i = 1:N
                    u = u - this.post_data.b_i_ell(i, ell) * sum(this.post_data.g_vecs(:, i)) * this.post_data.u_i_ell{i}(:, ell);
                end
            end
            tmp1 = this.opt_prob_interface.Apply_Misfit_Hessian(u, this.u_opt, this.z_opt);
            z_tmp = this.opt_prob_interface.Apply_Solution_Operator_z_Jacobian_Transpose(tmp1, this.z_opt);

            for ell = 1:N
                z_tmp = z_tmp + (this.post_data.state_grad' * this.post_data.u_ell(:, ell)) * (this.post_data.W_z_inv_Z(:, ell) - this.post_data.W_z_inv_z_opt);
                for i = 1:N
                    z_tmp = z_tmp - this.post_data.b_i_ell(i, ell) * (this.post_data.state_grad' * this.post_data.u_i_ell{i}(:, ell)) * (this.post_data.W_z_inv_Z * this.post_data.g_vecs(:, i) - sum(this.post_data.g_vecs(:, i)) * this.post_data.W_z_inv_z_opt);
                end
            end

            z_tmp = (1 / this.post_data.alpha_d) * z_tmp;
            z_pert = -this.md_hessian_analysis.Apply_RS_Hessian_Inverse(z_tmp, this.z_opt);
            z_update_mean = this.z_opt + z_pert;
        end

    end

end
