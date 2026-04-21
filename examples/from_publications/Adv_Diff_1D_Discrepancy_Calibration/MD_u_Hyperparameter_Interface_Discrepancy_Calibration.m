%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef MD_u_Hyperparameter_Interface_Discrepancy_Calibration < MD_u_Hyperparameter_Interface

    properties
        x
    end

    methods (Access = public)

        function [nodes] = Load_Spatial_Node_Data(this)
            nodes = cell(1, 1);
            nodes{1} = this.x;
        end

        function this = MD_u_Hyperparameter_Interface_Discrepancy_Calibration(x, data_centering)
            this@MD_u_Hyperparameter_Interface(false, data_centering, false);
            this.x = x;
        end

    end

end
