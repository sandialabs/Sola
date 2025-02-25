classdef MD_Prior_Sampling < handle

    properties
        data_interface
        u_prior_interface
        z_prior_interface
        u_opt
        z_opt
    end

    methods

        function this = MD_Prior_Sampling(data_interface, u_prior_interface, z_prior_interface)
            this.data_interface = data_interface;
            this.u_prior_interface = u_prior_interface;
            this.z_prior_interface = z_prior_interface;
            this.u_opt = this.data_interface.u_opt;
            this.z_opt = this.data_interface.z_opt;
        end

        function [delta_samples] = Prior_Discrepancy_Samples_at_z_opt(this, num_samps)
            delta_samples = this.u_prior_interface.Sample_with_Covariance_W_u_Inverse(num_samps) + this.data_interface.data_shift;
        end

        function [delta_samples, delta_zopt_samples] = Prior_Discrepancy_Samples(this, z, num_samps)
            Z = z - this.z_opt;
            Mz_Z = this.z_prior_interface.Apply_M_z(Z);
            Sigma = Mz_Z' * this.z_prior_interface.Apply_W_z_Inverse(Mz_Z);
            p = size(Z, 2);
            R = chol(Sigma);

            delta_zopt_samples = zeros(size(this.data_interface.D, 1), num_samps);
            delta_samples = cell(num_samps, 1);
            for k = 1:num_samps
                u_vec = this.u_prior_interface.Sample_with_Covariance_W_u_Inverse(p + 1);
                delta_samples{k} = u_vec(:, 1:p) * R + u_vec(:, p + 1) + this.data_interface.data_shift;
                delta_zopt_samples(:, k) = u_vec(:, p + 1) + this.data_interface.data_shift;
            end
        end

        function [delta_samples_z_opt, delta_samples_z_pert, z_pert] = Prior_Discrepancy_Samples_for_Visualization(this, num_samps, num_perts_init)
            E_z_inv_gsvd = E_z_Inv_GSVD(this.z_prior_interface, this.z_opt);
            [z_pert, ~, evals] = E_z_inv_gsvd.Compute_GSVD(num_perts_init, 10, 10);
            while evals(end) > .1
                num_perts_init = 2 * num_perts_init;
                [z_pert, ~, evals] = E_z_inv_gsvd.Compute_GSVD(num_perts_init, 10, 10);
            end
            I = find(evals>.1);
            z_pert = z_pert(:,I);
            evals = evals(I);
            num_perts = size(z_pert,2);
            delta_samples_z_opt = this.u_prior_interface.Sample_with_Covariance_W_u_Inverse(num_samps);

            scaling = .3 * sqrt(this.z_opt' * this.z_prior_interface.Apply_M_z(this.z_opt));
            z_pert = scaling * z_pert;

            delta_samples_z_pert = cell(num_perts, 1);
            for k = 1:num_perts
                delta_samples_z_pert{k} = scaling * sqrt(this.z_prior_interface.alpha_z) * evals(k) * this.u_prior_interface.Sample_with_Covariance_W_u_Inverse(num_samps);
            end
        end

    end

end
