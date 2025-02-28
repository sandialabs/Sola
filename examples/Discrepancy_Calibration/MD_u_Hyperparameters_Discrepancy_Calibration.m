classdef MD_u_Hyperparameters_Discrepancy_Calibration < MD_u_Hyperparameters

    properties
        x
    end

    methods (Access = public)

        function [nodes] = Load_Spatial_Node_Data(this)
            nodes = this.x;
        end

        function this = MD_u_Hyperparameters_Discrepancy_Calibration(data_interface, x, data_centering)
            this@MD_u_Hyperparameters(data_interface, false, data_centering);
            this.x = x;
        end

    end

end
