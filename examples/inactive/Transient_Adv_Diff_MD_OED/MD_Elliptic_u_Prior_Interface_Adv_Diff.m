classdef MD_Elliptic_u_Prior_Interface_Adv_Diff < MD_Elliptic_u_Prior_Interface

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

        function [u_out] = Apply_M_u(this, u_in)
            u_out = this.M * u_in;
        end

        function this = MD_Elliptic_u_Prior_Interface_Adv_Diff(alpha_u, sabl_opt)
            this@MD_Elliptic_u_Prior_Interface(alpha_u);

            S = sabl_opt.con.S;
            this.M = sabl_opt.con.M;

            this.E_u = (2.e-2) * S + this.M;

            num_sing_vals = sabl_opt.con.n_y;
            oversampling = 0;
            num_subspace_iters = 1;
            u_vec = zeros(sabl_opt.obj.n_y, 1);
            this.Compute_E_u_Inverse_GSVD(num_sing_vals, oversampling, num_subspace_iters, u_vec);
        end

    end

end
