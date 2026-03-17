classdef MD_u_Hyperparameter_Interface_Diff < MD_u_Hyperparameter_Interface

    properties
        x
    end

    methods (Access = public)

        function [nodes] = Load_Spatial_Node_Data(this)
            nodes = cell(1, 1);
            nodes{1} = this.x;
        end

        function this = MD_u_Hyperparameter_Interface_Diff(x)
            this@MD_u_Hyperparameter_Interface(false, false);
            this.x = x;
        end

    end

end
