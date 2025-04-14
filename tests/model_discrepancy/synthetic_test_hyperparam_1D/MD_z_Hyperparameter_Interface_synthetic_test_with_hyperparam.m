classdef MD_z_Hyperparameter_Interface_synthetic_test_with_hyperparam < MD_z_Hyperparameter_Interface

    properties
        x
    end

    methods (Access = public)

        function [nodes] = Load_Spatial_Node_Data(this)
            nodes = this.x;
        end

        function this = MD_z_Hyperparameter_Interface_synthetic_test_with_hyperparam(m)
            this@MD_z_Hyperparameter_Interface('spatial field');
            this.x = linspace(0, 1, m)';
        end

    end

end
