%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef MD_z_Hyperparameter_Interface_Discrepancy_Calibration < MD_z_Hyperparameter_Interface

    properties
        x
        con_lofi
    end

    methods (Access = public)

        function [nodes] = Load_Spatial_Node_Data(this)
            nodes = this.x;
        end

        function [u] = State_Solve(this, z)
            u = this.con_lofi.State_Solve(z);
        end

        function this = MD_z_Hyperparameter_Interface_Discrepancy_Calibration(num_state_solves, x, con_lofi)
            this@MD_z_Hyperparameter_Interface('spatial field', num_state_solves);
            this.x = x;
            this.con_lofi = con_lofi;
        end

    end

end
