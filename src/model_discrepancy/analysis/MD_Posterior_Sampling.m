classdef MD_Posterior_Sampling < handle

    properties
        data_interface
        u_prior_interface
        z_prior_interface
        post_data
        z_opt
    end

    methods

        function this = MD_Posterior_Sampling(data_interface, u_prior_interface, z_prior_interface)
            this.data_interface = data_interface;
            this.u_prior_interface = u_prior_interface;
            this.z_prior_interface = z_prior_interface;
            this.post_data = MD_Posterior_Data();
            this.z_opt = data_interface.z_opt;
        end

        function [] = Compute_Posterior_Data(this, alpha_d, num_samples)
            this.post_data.Compute_Posterior_Data(this.data_interface, this.u_prior_interface, this.z_prior_interface, alpha_d, this.z_opt, num_samples);
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
                    disp('Error in Posterior Discrepancy Samples: delta breve coeff < 0');
                end
                breve_coeff = sqrt(abs(tmp));
                delta_samples_k = delta_samples_k + breve_coeff * this.post_data.u_breve;

                delta_mean{k} = delta_mean_k;
                delta_samples{k} = delta_samples_k + delta_mean_k;

            end
        end

    end

end
