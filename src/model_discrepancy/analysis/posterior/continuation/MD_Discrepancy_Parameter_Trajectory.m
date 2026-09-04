%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef MD_Discrepancy_Parameter_Trajectory < Auxillary_Parameter_Trajectory

    properties
        sample_idx
    end

    methods

        function this = MD_Discrepancy_Parameter_Trajectory(num_continuation_steps, sample_index)
            this@Auxillary_Parameter_Trajectory(num_continuation_steps);
            this.sample_idx = sample_index;
        end

        function [t] = Get_Time(this, time_index)
            t = time_index / this.Get_Number_of_Timesteps();
        end

        function [sample_index] = Get_Sample_Index(this)
            sample_index = this.sample_idx;
        end

    end
end
