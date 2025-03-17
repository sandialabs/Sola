classdef MD_u_Hyperparameters_synthetic_test_with_hyperparam_auto < MD_u_Hyperparameters

    properties
        x
    end

    methods (Access = public)

        function [nodes] = Load_Spatial_Node_Data(this)
            nodes = cell(1,1);
            nodes{1} = this.x;
        end

        function this = MD_u_Hyperparameters_synthetic_test_with_hyperparam_auto(data_interface, m)
            this@MD_u_Hyperparameters(data_interface, false);
            this.x = linspace(0, 1, m)';
        end

    end

end
