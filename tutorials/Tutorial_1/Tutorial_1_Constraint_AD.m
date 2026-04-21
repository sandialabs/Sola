%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Tutorial_1_Constraint_AD < Constraint_AD

    methods (Access = public)

        function [c] = c_AD(this, u, z)
            c = [u(1) + u(2) - z(1)
                 z(1) * u(2) - z(2)
                 u(3)^3 - z(2)^2];
        end

    end
end
