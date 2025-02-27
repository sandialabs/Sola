classdef MD_z_Hyperparameters_multi_state_synthetic_test < MD_z_Hyperparameters

    properties
        x
    end

    methods (Access = public)

        function [nodes] = Load_Spatial_Node_Data(this)
            nodes = this.x;
        end

        function this = MD_z_Hyperparameters_multi_state_synthetic_test(data_interface, u_prior_interface, m)
            this@MD_z_Hyperparameters(data_interface, u_prior_interface);
            this.x = linspace(0, 1, m)';
        end

    end

end
