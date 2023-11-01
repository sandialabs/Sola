classdef HDSA_Sabl_MD_Interface_Elliptic_Prior_PDE_Test_Prob < HDSA_Sabl_MD_Interface_Elliptic_Prior

    properties
        E_u
        E_z
        E_d
        M
    end

    methods (Access = public)

        function [u_out] = Apply_E_u_Inverse(this, u_in)
            u_out = linsolve(this.E_u, u_in);
        end

        function [u_out] = Apply_E_u_Inverse_Transpose(this, u_in)
            u_out = linsolve(this.E_u', u_in);
        end

        function [u_out] = Apply_M_u(this, u_in)
            u_out = this.M * u_in;
        end

        function [u_out] = Apply_M_u_Inverse(this, u_in)
            u_out = linsolve(this.M, u_in);
        end

        function [z_out] = Apply_E_z_Inverse(this, z_in)
            z_out = linsolve(this.E_z, z_in);
        end

        function [z_out] = Apply_E_z_Inverse_Transpose(this, z_in)
            z_out = linsolve(this.E_z', z_in);
        end

        function [z_out] = Apply_M_z(this, z_in)
            z_out = this.M * z_in;
        end

        function [u_out] = Apply_E_d(this, u_in)
            u_out = this.E_d * u_in;
        end

        function [u_out] = Apply_E_d_Transpose(this, u_in)
            u_out = this.E_d' * u_in;
        end

        function [u_opt] = Load_Optimal_u(this)
            u_opt = load('u_opt.mat').u_opt;
        end

        function [z_opt] = Load_Optimal_z(this)
            z_opt = load('z_opt.mat').z_opt;
        end

        function [Z] = Load_Z_Data(this)
            Z = load('Z.mat').Z;
        end

        function [D] = Load_d_Data(this)
            D = load('D.mat').D;
        end

    end

    methods

        function this = HDSA_Sabl_MD_Interface_Elliptic_Prior_PDE_Test_Prob(con_opt_obj, alpha_u, alpha_z)
            this@HDSA_Sabl_MD_Interface_Elliptic_Prior(con_opt_obj, alpha_u, alpha_z);

            this.E_d = (1.e-6) * con_opt_obj.con.S + con_opt_obj.con.M;
            this.E_u = (1 * 10^-2) * con_opt_obj.con.S + con_opt_obj.con.M;
            this.E_z = (10^-3) * con_opt_obj.con.S + con_opt_obj.con.M;
            this.M = con_opt_obj.con.M;

            num_sing_vals = 200;
            oversampling = 0;
            num_subspace_iters = 1;
            u_vec = zeros(con_opt_obj.obj.m, 1);
            this.Compute_Elliptic_GSVD(num_sing_vals, oversampling, num_subspace_iters, u_vec);

        end

    end

end
