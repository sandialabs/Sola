classdef MD_u_Hyperparameter_Interface_Diff < MD_u_Hyperparameter_Interface

    properties
        x
        y
    end

    methods (Access = public)

        function [nodes] = Load_Spatial_Node_Data(this)
            nodes = cell(1, 1);
            nodes{1} = [this.x, this.y];
        end

        function this = MD_u_Hyperparameter_Interface_Diff(x, y, data_centering)
            this@MD_u_Hyperparameter_Interface(false, data_centering);
            this.x = x;
            this.y = y;
        end

    end

end
