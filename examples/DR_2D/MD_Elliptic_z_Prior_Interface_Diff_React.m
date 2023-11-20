classdef MD_Elliptic_z_Prior_Interface_Diff_React < MD_Elliptic_z_Prior_Interface

    properties
        E_z
        M
    end

    methods

        function [z_out] = Apply_E_z_Inverse(this, z_in)
            z_out = this.E_z \ z_in;
        end

        function [z_out] = Apply_E_z_Inverse_Transpose(this, z_in)
            z_out = this.E_z' \ z_in;
        end

        function [z_out] = Apply_M_z(this, z_in)
            z_out = this.M * z_in;
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
            z_out = this.M \ z_in;
        end

        function this = MD_Elliptic_z_Prior_Interface_Diff_React(alpha_z, sabl_opt)
            this@MD_Elliptic_z_Prior_Interface(alpha_z);

            S = sabl_opt.con.diff_react_lofi.pde_meshing.S;
            this.M = sabl_opt.con.diff_react_lofi.pde_meshing.M;
            this.E_z = (1.e-1) * S + this.M;
        end

    end

end
