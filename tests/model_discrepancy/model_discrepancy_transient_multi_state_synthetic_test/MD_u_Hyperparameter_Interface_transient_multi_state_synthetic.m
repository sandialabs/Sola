classdef MD_u_Hyperparameter_Interface_transient_multi_state_synthetic < MD_u_Hyperparameter_Interface

    properties
        n_y
        n_t
    end

    methods (Access = public)

        function [nodes] = Load_Spatial_Node_Data(this)
            nodes = cell(2, 1);
            nodes{1} = linspace(0, 1, this.n_y)';
            nodes{2} = linspace(0, 1, this.n_y)';
        end

        function [time_nodes] = Load_Time_Node_Data(this)
            time_nodes = linspace(0, 1, this.n_t)';
        end

        function this = MD_u_Hyperparameter_Interface_transient_multi_state_synthetic(n_y, n_t, component_id)
            this@MD_u_Hyperparameter_Interface(true, false, true, component_id);
            this.n_y = n_y;
            this.n_t = n_t;
        end

    end

end
