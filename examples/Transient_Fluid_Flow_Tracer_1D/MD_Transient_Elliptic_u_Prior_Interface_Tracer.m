classdef MD_Transient_Elliptic_u_Prior_Interface_Tracer < MD_Transient_Elliptic_u_Prior_Interface

    properties
        E_u
        M
    end

    methods

        function [u_out] = Apply_E_u_Inverse(this, u_in)
            u_out = linsolve(this.E_u, u_in);
        end

        function [u_out] = Apply_E_u_Inverse_Transpose(this, u_in)
            u_out = linsolve(this.E_u', u_in);
        end

        function [u_out] = Apply_Spatial_M_u(this, u_in)
            u_out = this.M * u_in;
        end

        function this = MD_Transient_Elliptic_u_Prior_Interface_Tracer(alpha_u, transient_prior_cov, sabl_opt)
            this@MD_Transient_Elliptic_u_Prior_Interface(alpha_u, transient_prior_cov);

            S = sabl_opt.con.S_z;
            this.M = sabl_opt.con.M_z;

            this.E_u = (2.e-2) * S + this.M;

            num_sing_vals = transient_prior_cov.n_y;
            oversampling = 0;
            num_subspace_iters = 1;
            u_vec = zeros(transient_prior_cov.n_y, 1);
            this.Compute_E_u_Inverse_GSVD(num_sing_vals, oversampling, num_subspace_iters, u_vec);
        end

    end

end
