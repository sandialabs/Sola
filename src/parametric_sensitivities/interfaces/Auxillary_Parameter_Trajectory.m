classdef Auxillary_Parameter_Trajectory < handle

    properties
        N
    end

    methods

        function this = Auxillary_Parameter_Trajectory(N)
            this.N = N;
        end

        function [N] = Get_Number_of_Timesteps(this)
            N = this.N;
        end

    end
end
