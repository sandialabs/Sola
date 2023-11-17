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

        function [z_samples] = Prior_z_Samples(this, num_samps)
            n = length(this.z_opt);
            Omega = randn(n, num_samps);
            z_samples = this.z_prior_interface.Apply_W_z_Inverse_Factor(Omega);
        end

        function [delta_samples] = Prior_Discrepancy_Samples_at_z_opt(this, num_samps)
            m = length(this.u_opt);
            Omega = randn(m, num_samps);
            delta_samples = this.u_prior_interface.Apply_W_u_Inverse_Factor(Omega);
        end

        function [delta_samples] = Prior_Discrepancy_Samples(this, z, num_samps)
            Z = z - this.z_opt;
            Sigma = Z' * this.z_prior_interface.Apply_W_z_Inverse(Z);
            m = length(this.u_opt);
            p = size(Z, 2);
            R = chol(Sigma);

            delta_samples = cell(num_samps, 1);
            for k = 1:num_samps
                u_vec = randn(m, p) * R + randn(m, 1);
                delta_samples{k} = this.u_prior_interface.Apply_W_u_Inverse_Factor(u_vec);
            end
        end

    end

end
