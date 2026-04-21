%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef MD_Elliptic_z_Prior_Interface_Transient_Test_Problem < MD_Elliptic_z_Prior_Interface

    properties
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

        function [z_out] = Apply_M_z(this, z_in)
            z_out = this.M * z_in;
        end

    end

    methods

        function this = MD_Elliptic_z_Prior_Interface_Transient_Test_Problem(alpha_z, sola_opt)
            this@MD_Elliptic_z_Prior_Interface(alpha_z);
            S = sola_opt.con.S;
            this.M = sola_opt.con.M;
            this.E_z = (3.e-2) * S + this.M;
        end

    end

end
