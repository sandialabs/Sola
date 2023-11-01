classdef Adv_Diff_HDSA < HDSA_Sabl_MD_Interface_Elliptic_Prior

    properties
        E_u
        E_z
        E_d
        M
    end

    methods

        function this = Adv_Diff_HDSA(con_opt_obj, alpha_u, alpha_z)
            this@HDSA_Sabl_MD_Interface_Elliptic_Prior(con_opt_obj, alpha_u, alpha_z);

            S = con_opt_obj.con.adv_diff.pde_meshing.S;
            this.M = con_opt_obj.con.adv_diff.pde_meshing.M;

            this.E_u = (5.e-1) * S + this.M;

            this.E_z = this.con_opt_obj.con.control_basis' * this.M * this.con_opt_obj.con.control_basis;

            this.E_d = (1.e-8) * S + this.M;

            num_sing_vals = 1000;
            oversampling = 0;
            num_subspace_iters = 1;
            u_vec = zeros(con_opt_obj.obj.m, 1);
            this.Compute_Elliptic_GSVD(num_sing_vals, oversampling, num_subspace_iters, u_vec);
        end

        function [u_out] = Apply_E_u_Inverse(this, u_in)
            u_out = this.E_u \ u_in;
        end

        function [u_out] = Apply_E_u_Inverse_Transpose(this, u_in)
            u_out = this.E_u' \ u_in;
        end

        function [u_out] = Apply_M_u(this, u_in)
            u_out = this.M * u_in;
        end

        function [u_out] = Apply_M_u_Inverse(this, u_in)
            u_out = this.M \ u_in;
        end

        function [z_out] = Apply_E_z_Inverse(this, z_in)
            z_out = this.E_z \ z_in;
        end

        function [z_out] = Apply_E_z_Inverse_Transpose(this, z_in)
            z_out = this.E_z' \ z_in;
        end

        function [z_out] = Apply_M_z(this, z_in)
            z_out = this.con_opt_obj.con.control_basis' * this.M * this.con_opt_obj.con.control_basis * z_in;
        end

        % This function must be implemented to enable Hessian GEVP
        function [z_out] = Apply_E_z(this, z_in)
            z_out = this.E_z * z_in;
        end

        % This function must be implemented to enable Hessian GEVP
        function [z_out] = Apply_E_z_Transpose(this, z_in)
            z_out = this.E_z' * z_in;
        end

        % This function must be implemented to enable Hessian GEVP
        function [z_out] = Apply_M_z_Inverse(this, z_in)
            z_out = (this.con_opt_obj.con.control_basis' * this.M * this.con_opt_obj.con.control_basis) \ z_in;
        end

        function [u_out] = Apply_E_d(this, u_in)
            u_out = this.E_d * u_in;
        end

        function [u_out] = Apply_E_d_Transpose(this, u_in)
            u_out = this.E_d' * u_in;
        end

        function [u_opt] = Load_Optimal_u(this)
            u_opt = load('Optimization_Results.mat').u_lofi;
        end

        function [z_opt] = Load_Optimal_z(this)
            z_opt = load('Optimization_Results.mat').z_lofi;
        end

        function [Z] = Load_Z_Data(this)
            Z = load('Optimization_Results.mat').Z;
        end

        function [D] = Load_d_Data(this)
            D = load('Optimization_Results.mat').D;
        end

    end

end
