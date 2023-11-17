classdef MD_Elliptic_u_Prior_Interface_Transient_Test_Problem < MD_Elliptic_u_Prior_Interface

    properties
        E_u
        E_d
        M
    end

    methods (Access = public)

        function [u_out] = Apply_E_u_Inverse(this, u_in)
            u_out = linsolve(this.E_u, u_in);
        end

        function [u_out] = Apply_E_u_Inverse_Transpose(this, u_in)
            u_out = linsolve(this.E_u, u_in);
        end

        function [u_out] = Apply_M_u(this, u_in)
            u_out = this.M * u_in;
        end

        function [u_out] = Apply_M_u_Inverse(this, u_in)
            u_out = linsolve(this.M, u_in);
        end

        function [u_out] = Apply_E_d(this, u_in)
            u_out = this.E_d * u_in;
        end

        function [u_out] = Apply_E_d_Transpose(this, u_in)
            u_out = this.E_d' * u_in;
        end

    end

    methods

        function this = MD_Elliptic_u_Prior_Interface_Transient_Test_Problem(alpha_u, transient_prior_cov, sabl_opt)
            this@MD_Elliptic_u_Prior_Interface(alpha_u, transient_prior_cov);

            S = sabl_opt.con.S;
            this.M = sabl_opt.con.M;

            this.E_u = (2.e-2) * S + this.M;

            this.E_d = (1.e-8) * S + this.M;

            num_sing_vals = sabl_opt.con.n_y;
            oversampling = 0;
            num_subspace_iters = 1;
            u_vec = zeros(sabl_opt.con.n_y, 1);
            this.Compute_E_u_Inverse_GSVD(num_sing_vals, oversampling, num_subspace_iters, u_vec);
        end

    end

end
