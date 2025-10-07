classdef MD_u_Hyperparameter_Interface_synthetic_test_lumped_mass < MD_u_Hyperparameter_Interface

    properties
        x
    end

    methods (Access = public)

        function [nodes] = Load_Spatial_Node_Data(this)
            nodes = cell(1, 1);
            nodes{1} = this.x;
        end

        function this = MD_u_Hyperparameter_Interface_synthetic_test_lumped_mass(m)
            this@MD_u_Hyperparameter_Interface(false, true);
            this.x = linspace(0, 1, m)';
        end

    end

end
