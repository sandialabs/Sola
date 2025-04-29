classdef MD_u_Hyperparameter_Interface_synthetic_test_transient < MD_u_Hyperparameter_Interface

    properties
        x
        t
    end

    methods (Access = public)

        function [nodes] = Load_Spatial_Node_Data(this)
            nodes = cell(1, 1);
            nodes{1} = this.x;
        end

        function [nodes] = Load_Time_Node_Data(this)
            nodes = this.t;
        end

        function this = MD_u_Hyperparameter_Interface_synthetic_test_transient(n_y,n_t,T)
            this@MD_u_Hyperparameter_Interface(true);
            this.x = linspace(0, 1, n_y)';
            this.t = linspace(0,T,n_t)';
        end

    end

end
