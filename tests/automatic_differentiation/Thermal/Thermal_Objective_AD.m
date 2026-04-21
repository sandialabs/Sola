%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Thermal_Objective_AD < Objective_AD

    properties
        con
    end

    methods (Access = public)

        function [val] = J_AD(this, u, z)
            val = (1 / 2) * u' * this.con.M * u + (1 / 2) * z' * this.con.S * z;
        end

    end

    methods (Access = public)

        function this = Thermal_Objective_AD(n_u, n_z, con)
            this@Objective_AD(n_u, n_z);
            this.con = con;
        end

    end
end
