classdef MD_Transient_Elliptic_u_Prior_Interface_Transient_ADR_2D < MD_Transient_Elliptic_u_Prior_Interface

    properties
        M_u
        E_u
    end

    methods (Access = public)

        function [u_out] = Apply_E_u_Inverse(this, u_in)
            u_out = this.E_u \ u_in;
        end

        function [u_out] = Apply_E_u_Inverse_Transpose(this, u_in)
            u_out = (this.E_u') \ u_in;
        end

        function [u_out] = Apply_Spatial_M_u(this, u_in)
            u_out = this.M_u * u_in;
        end

        function this = MD_Transient_Elliptic_u_Prior_Interface_Transient_ADR_2D(alpha_u, transient_prior_cov, M, S)
            this@MD_Transient_Elliptic_u_Prior_Interface(alpha_u, transient_prior_cov);
            this.M_u = M;
            this.E_u = (1.e-3) * S + M; % Need to look at this more closely

            num_sing_vals = 100; % Need to look at this more closely
            oversampling = 20;
            num_subspace_iters = 1;
            u_vec = zeros(size(M, 1), 1);
            this.Compute_E_u_Inverse_GSVD(num_sing_vals, oversampling, num_subspace_iters, u_vec);

        end

    end

end
