classdef MD_z_Prior_Interface_Transient_ADR_2D < MD_z_Prior_Interface

    properties
        W_z
        R
    end

    methods (Access = public)

        function [z_out] = Apply_W_z_Inverse(this, z_in)
            z_out = linsolve(this.W_z, z_in);
        end

        function [z_out] = Sample_with_Covariance_W_z_Inverse(this, num_samples)
            omega = randn(size(this.R, 1), num_samples);
            z_out = linsolve(this.R, omega);
        end

        function [z_out] = Apply_W_z(this, z_in)
            z_out = this.W_z * z_in;
        end

        function this = MD_z_Prior_Interface_Transient_ADR_2D(transient_prior_cov, n_q)
            this.W_z = kron(transient_prior_cov.E_t, eye(n_q));
            this.R = chol(this.W_z);
        end

    end

end
