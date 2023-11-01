classdef HDSA_MD_Prior_Sampling < handle

    properties
        md_interface
        u_opt
        z_opt
    end

    methods

        function this = HDSA_MD_Prior_Sampling(md_interface)
            this.md_interface = md_interface;
            this.u_opt = this.md_interface.Load_Optimal_u();
            this.z_opt = this.md_interface.Load_Optimal_z();
        end

        function [z_samples] = Prior_z_Samples(this, num_samps)
            n = length(this.z_opt);
            Omega = randn(n, num_samps);
            z_samples = this.md_interface.Apply_W_z_Inverse_Factor(Omega);
        end

        function [delta_samples] = Prior_Discrepancy_Samples_at_z_opt(this, num_samps)
            m = length(this.u_opt);
            Omega = randn(m, num_samps);
            delta_samples = this.md_interface.Apply_W_u_Inverse_Factor(Omega);
        end

        function [delta_samples] = Prior_Discrepancy_Samples(this, z, num_samps)
            Z = z - this.z_opt;
            Sigma = Z' * this.md_interface.Apply_W_z_Inverse(Z);
            m = length(this.u_opt);
            p = size(Z, 2);
            R = chol(Sigma);

            delta_samples = cell(num_samps, 1);
            for k = 1:num_samps
                u_vec = randn(m, p) * R + randn(m, 1);
                delta_samples{k} = this.md_interface.Apply_W_u_Inverse_Factor(u_vec);
            end
        end

    end

end
