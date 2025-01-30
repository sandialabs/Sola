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

        function [delta_mean, delta_samples] = Posterior_Discrepancy_Samples(this, z)
            m = size(this.post_data.u_ell, 1);
            p = size(z, 2);
            delta_mean = cell(p, 1);
            delta_samples = cell(p, 1);

            for k = 1:p
                dz = z(:, k) - this.z_opt;

                M_z_dz = this.z_prior_interface.Apply_M_z(dz);
                coeffs = 1 + (this.post_data.W_z_inv_M_z_Z - this.post_data.W_z_inv_M_z_z_opt)' * M_z_dz;
                delta_mean_k = this.post_data.u_ell * coeffs;

                delta_samples_k = zeros(m, this.post_data.num_samples);
                for i = 1:this.post_data.N
                    sgi = sum(this.post_data.g_vecs(:, i));
                    W_z_inv_M_z_yi = this.post_data.W_z_inv_M_z_Z * this.post_data.g_vecs(:, i) - sgi * this.post_data.W_z_inv_M_z_z_opt;
                    coeffs = this.post_data.b_i_ell(i, :)' * (sgi + W_z_inv_M_z_yi' * M_z_dz);
                    delta_mean_k = delta_mean_k - this.post_data.u_i_ell{i} * coeffs;

                    coeff = (1 / sqrt(this.post_data.Mu(i, i))) * (sgi + W_z_inv_M_z_yi' * M_z_dz);
                    delta_samples_k = delta_samples_k + coeff * this.post_data.ui_hat{i};
                end
                delta_mean_k = (1 / this.post_data.alpha_d) * delta_mean_k;
                delta_samples_k = sqrt(this.post_data.alpha_d) * delta_samples_k;

                W_z_inv_M_z_dz = this.z_prior_interface.Apply_W_z_Inverse(M_z_dz);
                tmp = M_z_dz' * W_z_inv_M_z_dz - W_z_inv_M_z_dz' * this.post_data.M_z_Zc * linsolve(this.post_data.Zc_M_z_W_z_inv_M_z_Zc, this.post_data.M_z_Zc' * W_z_inv_M_z_dz);
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
