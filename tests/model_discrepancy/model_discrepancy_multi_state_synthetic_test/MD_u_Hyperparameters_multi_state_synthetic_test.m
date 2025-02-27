classdef MD_u_Hyperparameters_multi_state_synthetic_test < MD_u_Hyperparameters

    properties
        x
    end

    methods (Access = public)

        function [nodes] = Load_Spatial_Node_Data(this)
            nodes = this.x;
        end

        function this = MD_u_Hyperparameters_multi_state_synthetic_test(data_interface, component_id, m)
            this@MD_u_Hyperparameters(data_interface, false, [], [], component_id);
            this.x = linspace(0, 1, m)';
        end

    end

end
