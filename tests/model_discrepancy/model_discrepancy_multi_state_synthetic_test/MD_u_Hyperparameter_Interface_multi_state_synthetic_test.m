classdef MD_u_Hyperparameter_Interface_multi_state_synthetic_test < MD_u_Hyperparameter_Interface

    properties
        x
    end

    methods (Access = public)

        function [nodes] = Load_Spatial_Node_Data(this)
            nodes = cell(2,1);
            nodes{1} = this.x;
            nodes{2} = this.x;
        end

        function this = MD_u_Hyperparameter_Interface_multi_state_synthetic_test(component_id, m)
            this@MD_u_Hyperparameter_Interface(false, [], [], component_id);
            this.x = linspace(0, 1, m)';
        end

    end

end
