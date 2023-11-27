classdef MD_Elliptic_z_Prior_Interface_PDE_Test_Problem < MD_Elliptic_z_Prior_Interface

    properties
        m
        x
        M
        E_z
    end

    methods (Access = public)

        function [z_out] = Apply_E_z_Inverse(this, z_in)
            z_out = linsolve(this.E_z, z_in);
        end

        function [z_out] = Apply_E_z_Inverse_Transpose(this, z_in)
            z_out = linsolve(this.E_z', z_in);
        end

        function [z_out] = Apply_E_z(this, z_in)
            z_out = this.E_z * z_in;
        end

        function [z_out] = Apply_E_z_Transpose(this, z_in)
            z_out = this.E_z' * z_in;
        end

        function [z_out] = Apply_M_z(this, z_in)
            z_out = this.M * z_in;
        end

        function [z_out] = Apply_M_z_Inverse(this, z_in)
            z_out = linsolve(this.M, z_in);
        end

    end

    methods

        function this = MD_Elliptic_z_Prior_Interface_PDE_Test_Problem(alpha_z, sabl_opt)
            this@MD_Elliptic_z_Prior_Interface(alpha_z);
            this.E_z = (10^-3) * sabl_opt.con.S + sabl_opt.con.M;
            this.M = sabl_opt.con.M;
        end

    end

end
