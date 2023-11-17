classdef MD_Elliptic_u_Prior_Interface_Adv_Diff < MD_Elliptic_u_Prior_Interface

    properties
        E_u
        M_u
        E_d
    end

    methods (Access = public)

        function [u_out] = Apply_E_u_Inverse(this, u_in)
            u_out = linsolve(this.E_u, u_in);
        end

        function [u_out] = Apply_E_u_Inverse_Transpose(this, u_in)
            u_out = linsolve(this.E_u, u_in);
        end

        function [u_out] = Apply_M_u(this, u_in)
            u_out = this.M_u * u_in;
        end

        function [u_out] = Apply_M_u_Inverse(this, u_in)
            u_out = linsolve(this.M_u, u_in);
        end

        function [u_out] = Apply_E_d(this, u_in)
            u_out = this.E_d * u_in;
        end

        function [u_out] = Apply_E_d_Transpose(this, u_in)
            u_out = this.E_d' * u_in;
        end

    end

    methods

        function this = MD_Elliptic_u_Prior_Interface_Adv_Diff(alpha_u, sabl_opt)
            this@MD_Elliptic_u_Prior_Interface(alpha_u);

            this.E_d = (1.e-6) * sabl_opt.con.S + sabl_opt.con.M;
            this.E_u = (1 * 10^-2) * sabl_opt.con.S + sabl_opt.con.M;
            this.M_u = sabl_opt.con.M;

            num_sing_vals = 200;
            oversampling = 0;
            num_subspace_iters = 1;
            u_vec = zeros(sabl_opt.obj.m, 1);
            this.Compute_E_u_Inverse_GSVD(num_sing_vals, oversampling, num_subspace_iters, u_vec);

        end

    end

end
