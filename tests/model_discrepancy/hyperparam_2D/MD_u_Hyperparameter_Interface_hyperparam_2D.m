%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef MD_u_Hyperparameter_Interface_hyperparam_2D < MD_u_Hyperparameter_Interface

    properties
        x
        y
    end

    methods (Access = public)

        function [spatial_nodes] = Load_Spatial_Node_Data(this)
            spatial_nodes = cell(1, 1);
            spatial_nodes{1} = [this.x, this.y];
        end

        function this = MD_u_Hyperparameter_Interface_hyperparam_2D(x, y)
            this@MD_u_Hyperparameter_Interface(false);
            this.x = x;
            this.y = y;
        end

    end

end
