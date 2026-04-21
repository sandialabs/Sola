%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Auxillary_Parameter_Trajectory < handle

    properties
        N
    end

    methods

        function this = Auxillary_Parameter_Trajectory(N)
            arguments
                N (1, 1) {mustBeNumeric}
            end
            this.N = N;
        end

        function [N] = Get_Number_of_Timesteps(this)
            N = this.N;
        end

    end
end
