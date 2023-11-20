classdef  MD_Elliptic_z_Prior_Interface_Adv_Diff < MD_Elliptic_z_Prior_Interface

    properties
        E_z
        M
        control_basis
    end

    methods

        function [z_out] = Apply_E_z_Inverse(this, z_in)
            z_out = this.E_z \ z_in;
        end

        function [z_out] = Apply_E_z_Inverse_Transpose(this, z_in)
            z_out = this.E_z' \ z_in;
        end

        function [z_out] = Apply_M_z(this, z_in)
            z_out = this.control_basis' * this.M * this.control_basis * z_in;
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
            z_out = (this.control_basis' * this.M * this.control_basis) \ z_in;
        end

        function this = MD_Elliptic_z_Prior_Interface_Adv_Diff(alpha_z, sabl_opt)
            this@MD_Elliptic_z_Prior_Interface(alpha_z);

            this.M = sabl_opt.con.adv_diff.pde_meshing.M;
            this.E_z = sabl_opt.con.control_basis' * this.M * sabl_opt.con.control_basis;
            this.control_basis = sabl_opt.con.control_basis;
        end

    end

end
