%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Mean-to-sample trajectory extrapolated past the sample: the discrepancy at
% pseudo-time t in [0, 1] is  mean + scale * t * (sample - mean).  scale = 1 is
% the ordinary mean-to-sample trajectory; scale > 1 ends at the inflated
% "sample" mean + scale * (sample - mean), which is a controlled way to make
% the increment (and therefore the continuation) harder.  The sensitivity
% operators multiply the pseudo-time derivative by Get_Scale().

classdef MD_Discrepancy_Extrapolated_Trajectory < MD_Discrepancy_Mean_To_Sample_Trajectory

    properties
        scale
    end

    methods

        function this = MD_Discrepancy_Extrapolated_Trajectory(num_continuation_steps, sample_index, scale)
            this@MD_Discrepancy_Mean_To_Sample_Trajectory(num_continuation_steps, sample_index);
            this.scale = scale;
        end

        function [t] = Get_Time(this, time_index)
            t = this.scale * time_index / this.Get_Number_of_Timesteps();
        end

        function [scale] = Get_Scale(this)
            scale = this.scale;
        end

    end
end
