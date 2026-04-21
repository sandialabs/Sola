%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef MD_Transient_Elliptic_u_Prior_Interface_Transient_Test_Problem < MD_Transient_Elliptic_u_Prior_Interface

    properties
        E_u
        M
    end

    methods (Access = public)

        function [u_out] = Apply_E_u_Inverse(this, u_in)
            u_out = linsolve(this.E_u, u_in);
        end

        function [u_out] = Apply_E_u_Inverse_Transpose(this, u_in)
            u_out = linsolve(this.E_u, u_in);
        end

        function [u_out] = Apply_Spatial_M_u(this, u_in)
            u_out = this.M * u_in;
        end

    end

    methods

        function this = MD_Transient_Elliptic_u_Prior_Interface_Transient_Test_Problem(alpha_u, transient_prior_cov, sabl_opt, num_sing_vals)
            this@MD_Transient_Elliptic_u_Prior_Interface(alpha_u, transient_prior_cov);

            S = sabl_opt.con.S;
            this.M = sabl_opt.con.M;

            this.E_u = (2.e-2) * S + this.M;

            oversampling = sabl_opt.con.n_y - num_sing_vals;
            num_subspace_iters = 1;
            u_vec = zeros(sabl_opt.con.n_y, 1);
            this.Compute_E_u_Inverse_GSVD(num_sing_vals, oversampling, num_subspace_iters, u_vec);
        end

    end

end
