%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Thermal_Objective < Objective

    properties
        con
    end

    methods (Access = public)

        function [val, grad_u, grad_z] = J(this, u, z)
            val = (1 / 2) * u' * this.con.M * u + (1 / 2) * z' * this.con.S * z;
            grad_u = this.con.M * u;
            grad_z = this.con.S * z;
        end

        function [Mv] = J_uu_Apply(this, v, u, z)
            Mv = this.con.M * v;
        end

        function [Mv] = J_uz_Apply(this, v, u, z)
            Mv = 0 * u;
        end

        function [Mv] = J_zu_Apply(this, v, u, z)
            Mv = 0 * z;
        end

        function [Mv] = J_zz_Apply(this, v, u, z)
            Mv = this.con.S * v;
        end

    end

    methods (Access = public)

        function this = Thermal_Objective(con)
            this.con = con;
        end

    end
end
