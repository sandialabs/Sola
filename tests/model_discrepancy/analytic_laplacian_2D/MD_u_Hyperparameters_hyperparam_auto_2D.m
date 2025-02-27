classdef MD_u_Hyperparameters_hyperparam_auto_2D < MD_u_Hyperparameters

    properties
        x
        y
    end

    methods (Access = public)

        function [nodes] = Load_Spatial_Node_Data(this)
            nodes = [this.x, this.y];
        end

        function this = MD_u_Hyperparameters_hyperparam_auto_2D(data_interface, x, y)
            this@MD_u_Hyperparameters(data_interface, false);
            this.x = x;
            this.y = y;
        end

    end

end
