classdef MD_z_Hyperparameter_Interface_Transient_Test_Problem < MD_z_Hyperparameter_Interface

    properties
        x
    end

    methods (Access = public)

        function [nodes] = Load_Spatial_Node_Data(this)
            nodes = this.x;
        end

        function this = MD_z_Hyperparameter_Interface_Transient_Test_Problem(x)
            this@MD_z_Hyperparameter_Interface('spatial field');
            this.x = x;
        end

    end

end
