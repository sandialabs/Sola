%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Thermal_Objective_AD < Dynamic_Objective_AD

    properties

    end

    methods (Access = public)

        function [val] = g_AD(this, y, t)
            val = (1 / 2) * (y' * y);
        end

        function [val] = R_AD(this, z)
            val = (1 / 2) * (z' * z);
        end

    end

    methods (Access = public)

        function this = Thermal_Objective_AD(m, n, T, N)
            this@Dynamic_Objective_AD(m, n, T, N);
        end

    end

end
