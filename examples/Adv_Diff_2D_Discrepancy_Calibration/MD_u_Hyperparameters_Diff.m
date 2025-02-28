classdef MD_u_Hyperparameters_Diff < MD_u_Hyperparameters

    properties
        x
        y
    end

    methods (Access = public)

        function [nodes] = Load_Spatial_Node_Data(this)
            nodes = [this.x, this.y];
        end

        function this = MD_u_Hyperparameters_Diff(data_interface, x, y, data_centering)
            this@MD_u_Hyperparameters(data_interface, false, data_centering);
            this.x = x;
            this.y = y;
        end

    end

end
