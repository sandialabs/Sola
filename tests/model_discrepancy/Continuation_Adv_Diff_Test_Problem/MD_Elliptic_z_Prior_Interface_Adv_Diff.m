classdef MD_Elliptic_z_Prior_Interface_Adv_Diff< MD_Elliptic_z_Prior_Interface

    properties
        E_z
        M_z
    end

    methods (Access = public)

        function [z_out] = Apply_E_z_Inverse(this, z_in)
            z_out = linsolve(this.E_z, z_in);
        end

        function [z_out] = Apply_E_z_Inverse_Transpose(this, z_in)
            z_out = linsolve(this.E_z', z_in);
        end

        function [z_out] = Apply_M_z(this, z_in)
            z_out = this.M_z * z_in;
        end

    end

    methods

        function this = MD_Elliptic_z_Prior_Interface_Adv_Diff(alpha_z, sabl_opt)
            this@MD_Elliptic_z_Prior_Interface(alpha_z);
            this.E_z = (10^-3) * sabl_opt.con.S + sabl_opt.con.M;
            this.M_z = sabl_opt.con.M;
        end

    end

end
