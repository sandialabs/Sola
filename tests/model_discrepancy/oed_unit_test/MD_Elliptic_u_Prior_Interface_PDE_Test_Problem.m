classdef MD_Elliptic_u_Prior_Interface_PDE_Test_Problem < MD_Elliptic_u_Prior_Interface

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

        function [u_out] = Apply_M_u(this, u_in)
            u_out = this.M * u_in;
        end

    end

    methods

        function this = MD_Elliptic_u_Prior_Interface_PDE_Test_Problem(alpha_u, sabl_opt)
            this@MD_Elliptic_u_Prior_Interface(alpha_u);

            this.E_u = (1 * 10^-2) * sabl_opt.con.S + sabl_opt.con.M;
            this.M = sabl_opt.con.M;

            num_sing_vals = 200;
            oversampling = 0;
            num_subspace_iters = 1;
            u_vec = zeros(sabl_opt.con.m, 1);
            this.Compute_E_u_Inverse_GSVD(num_sing_vals, oversampling, num_subspace_iters, u_vec);

        end

    end

end
