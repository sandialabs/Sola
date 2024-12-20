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
            delta_samples = this.u_prior_interface.Sample_with_Covariance_W_u_Inverse(num_samps);
        end

        function [delta_samples] = Prior_Discrepancy_Samples(this, z, num_samps)
            Z = z - this.z_opt;
            Sigma = Z' * this.z_prior_interface.Apply_W_z_Inverse(Z);
            p = size(Z, 2);
            R = chol(Sigma);

            delta_samples = cell(num_samps, 1);
            for k = 1:num_samps
                u_vec = this.u_prior_interface.Sample_with_Covariance_W_u_Inverse(p + 1);
                delta_samples{k} = u_vec(:, 1:p) * R + u_vec(:, p + 1);
            end
        end

        function [sample_norms,discrepancy_data_norms] = Compute_Prior_Discrepancy_Norms(this,num_samples)
            z_opt_norm = sqrt(this.data_interface.z_opt'*this.z_prior_interface.Apply_M_z(this.data_interface.z_opt));
            z_samples = this.z_prior_interface.Sample_with_Covariance_W_z_Inverse(num_samples);
            for k = 1:num_samples
                z_samples(:,k) = z_opt_norm*z_samples(:,k)/sqrt(z_samples(:,k)'*this.z_prior_interface.Apply_M_z(z_samples(:,k)));
            end
            prior_delta_samples = this.Prior_Discrepancy_Samples(z_samples,num_samples);
            sample_norms = zeros(num_samples,num_samples);
            for i = 1:num_samples
                tmp = prior_delta_samples{i}'*this.u_prior_interface.Apply_M_u(prior_delta_samples{i});
                sample_norms(:,i) = sqrt(diag(tmp));
            end
            discrepancy_data_norms = sqrt(diag(this.data_interface.D'*this.u_prior_interface.Apply_M_u(this.data_interface.D)));
        end

    end

end
