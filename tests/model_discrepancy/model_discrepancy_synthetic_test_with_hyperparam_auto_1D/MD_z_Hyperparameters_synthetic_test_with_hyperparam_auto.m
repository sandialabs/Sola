classdef MD_z_Hyperparameters_synthetic_test_with_hyperparam_auto < MD_z_Hyperparameters

    properties
        x
    end

    methods (Access = public)

        function [nodes] = Load_Spatial_Node_Data(this)
            nodes = this.x;
        end

        function this = MD_z_Hyperparameters_synthetic_test_with_hyperparam_auto(data_interface, u_prior_interface, m)
            this@MD_z_Hyperparameters(data_interface, u_prior_interface, 'spatial field');
            this.x = linspace(0, 1, m)';
        end

    end

end
