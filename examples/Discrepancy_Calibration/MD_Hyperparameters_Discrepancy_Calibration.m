classdef MD_Hyperparameters_Discrepancy_Calibration < MD_Hyperparameters

    properties
        x
    end

    methods (Access = public)

        function [nodes] = Load_Node_Data(this)
            nodes = this.x;
        end

        function this = MD_Hyperparameters_Discrepancy_Calibration(data_interface,x)
            this@MD_Hyperparameters(data_interface);
            this.x = x;
        end

    end

end