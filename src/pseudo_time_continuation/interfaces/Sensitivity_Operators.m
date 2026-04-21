%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Sensitivity_Operators < handle

    properties

    end

    methods (Abstract, Access = public)

        [grad, val] = Gradient(this, z, theta_traj, time_index)

        [z_out] = Apply_Hessian(this, z_in, z, theta_traj, time_index)

        [z_out] = Apply_B(this, z, theta_traj, time_index)

    end

    methods

        function this = Sensitivity_Operators()

        end

    end
end
