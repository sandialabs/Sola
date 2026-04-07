classdef MD_Posterior_Sampling < handle

    properties
        data_interface
        u_prior_interface
        z_prior_interface
        post_data
        z_opt

        % NEW
        Mz_Wz_inv_Mz_Z_minus_z_opt
        Mz_Wz_inv_Mz_yi
        si
    end

    methods

        function this = MD_Posterior_Sampling(data_interface, u_prior_interface, z_prior_interface)
            arguments
                data_interface MD_Data_Interface
                u_prior_interface MD_u_Prior_Interface
                z_prior_interface MD_z_Prior_Interface
            end
            this.data_interface = data_interface;
            this.u_prior_interface = u_prior_interface;
            this.z_prior_interface = z_prior_interface;
            this.post_data = MD_Posterior_Data();
            this.z_opt = data_interface.z_opt;
        end

        function [] = Compute_Posterior_Data(this, alpha_d, num_samples)
            this.post_data.Compute_Posterior_Data(this.data_interface, this.u_prior_interface, this.z_prior_interface, alpha_d, num_samples);
            this.Mz_Wz_inv_Mz_Z_minus_z_opt = this.post_data.Mz_Wz_inv_Mz_Z - this.post_data.Mz_Wz_inv_Mz_z_opt;
            this.Mz_Wz_inv_Mz_yi = 0 * this.Mz_Wz_inv_Mz_Z_minus_z_opt;
            this.si = zeros(1, this.post_data.N);
            for i = 1:this.post_data.N
                this.Mz_Wz_inv_Mz_yi(:, i) = this.post_data.Mz_Wz_inv_Mz_Z * this.post_data.g_vecs(:, i) - ...
                    sum(this.post_data.g_vecs(:, i)) * this.post_data.Mz_Wz_inv_Mz_z_opt;
                this.si(i) = sum(this.post_data.g_vecs(:, i)) - this.z_opt' * this.Mz_Wz_inv_Mz_yi(:, i);
            end
        end

        function [delta_mean, delta_samples] = Posterior_Discrepancy_Samples(this, z)
            m = size(this.post_data.u_ell, 1);
            p = size(z, 2);
            delta_mean = cell(p, 1);
            delta_samples = cell(p, 1);

            for k = 1:p
                dz = z(:, k) - this.z_opt;

                Mz_dz = this.z_prior_interface.Apply_M_z(dz);
                coeffs = 1 + (this.post_data.Wz_inv_Mz_Z - this.post_data.Wz_inv_Mz_z_opt)' * Mz_dz;
                delta_mean_k = this.post_data.u_ell * coeffs;

                delta_samples_k = zeros(m, this.post_data.num_samples);
                for i = 1:this.post_data.N
                    sgi = sum(this.post_data.g_vecs(:, i));
                    Wz_inv_Mz_yi = this.post_data.Wz_inv_Mz_Z * this.post_data.g_vecs(:, i) - sgi * this.post_data.Wz_inv_Mz_z_opt;
                    coeffs = this.post_data.b_i_ell(i, :)' * (sgi + Wz_inv_Mz_yi' * Mz_dz);
                    delta_mean_k = delta_mean_k - this.post_data.u_i_ell{i} * coeffs;

                    coeff = (1 / sqrt(this.post_data.Mu(i, i))) * (sgi + Wz_inv_Mz_yi' * Mz_dz);
                    delta_samples_k = delta_samples_k + coeff * this.post_data.ui_hat{i};
                end
                delta_mean_k = (1 / this.post_data.alpha_d) * delta_mean_k + this.data_interface.data_shift;
                delta_samples_k = sqrt(this.post_data.alpha_d) * delta_samples_k;

                Wz_inv_Mz_dz = this.z_prior_interface.Apply_W_z_Inverse(Mz_dz);
                tmp = Mz_dz' * Wz_inv_Mz_dz - Wz_inv_Mz_dz' * this.post_data.Mz_Zc * linsolve(this.post_data.Zc_Mz_Wz_inv_Mz_Zc, this.post_data.Mz_Zc' * Wz_inv_Mz_dz);
                if tmp < -1.e-11
                    disp('Error in Posterior Discrepancy Samples: delta breve coeff < 0');
                end
                breve_coeff = sqrt(abs(tmp));
                delta_samples_k = delta_samples_k + breve_coeff * this.post_data.u_breve;

                delta_mean{k} = delta_mean_k;
                delta_samples{k} = delta_samples_k + delta_mean_k;

            end
        end

        % ------------------------------------------------------------
        % Discrepancy kernels: Mean
        % ------------------------------------------------------------

        function [u_out] = Discrepancy_Evaluation_Mean(this, z_n)
            N = this.post_data.N;
            u_out = zeros(size(this.data_interface.u_opt));
            for ell = 1:N
                coeff = this.post_data.a_ell(ell) + z_n' * this.Mz_Wz_inv_Mz_Z_minus_z_opt(:, ell);
                u_out = u_out + coeff * this.post_data.u_ell(:, ell);
                for i = 1:N
                    coeff = this.post_data.b_i_ell(i, ell) * (this.si(i) + z_n' * this.Mz_Wz_inv_Mz_yi(:, i));
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

        function [z_out] = Apply_Discrepancy_z_Jacobian_transpose_Mean(this, u_in)
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

        function [u_out] = Discrepancy_Evaluation_Sample(this, z_n, sample_idx)
            u_out_mean = this.Discrepancy_Evaluation_Mean(z_n);
            dz = z_n - this.z_opt;

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

        function [u_out] = Apply_Discrepancy_z_Jacobian_Sample(this, z_n, z_in, sample_idx)
            % Note: z_n is needed since sampling is nonlinear in z via gamma(z) in delta_breve
            u_out_mean = this.Apply_Discrepancy_z_Jacobian_Mean(z_in);

            u = zeros(size(u_out_mean));
            for i = 1:this.post_data.N
                coeff = (1 / sqrt(this.post_data.Mu(i, i))) * (this.Mz_Wz_inv_Mz_yi(:, i)' * z_in);
                u = u + coeff * this.post_data.ui_hat{i}(:, sample_idx);
            end
            u = sqrt(this.post_data.alpha_d) * u;

            Mz_dz = this.z_prior_interface.Apply_M_z(z_n - this.z_opt);
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

        function [z_out] = Apply_Discrepancy_z_Jacobian_transpose_Sample(this, z_n, u_in, sample_idx)
            % Note: z_n is needed since sampling is nonlinear in z via gamma(z) in delta_breve
            z_out_mean = this.Apply_Discrepancy_z_Jacobian_transpose_Mean(u_in);

            z = zeros(size(z_out_mean));
            for i = 1:this.post_data.N
                ui_hat_idx = this.post_data.ui_hat{i}(:, sample_idx);
                coeff = (1 / sqrt(this.post_data.Mu(i, i))) * (ui_hat_idx' * u_in);
                z = z + coeff * this.Mz_Wz_inv_Mz_yi(:, i);
            end
            z = sqrt(this.post_data.alpha_d) * z;

            Mz_dz = this.z_prior_interface.Apply_M_z(z_n - this.z_opt);
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
            z = z + (u_breve_idx' * u_in) * breve_coeff_grad;

            z_out = z_out_mean + z;
        end

    end

end
