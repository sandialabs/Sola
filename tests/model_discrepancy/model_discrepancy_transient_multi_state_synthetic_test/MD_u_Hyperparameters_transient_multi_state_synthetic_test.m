classdef MD_u_Hyperparameters_transient_multi_state_synthetic_test < MD_u_Hyperparameters

    properties
        n_y
        n_t
    end

    methods (Access = public)

        function [spatial_nodes] = Load_Spatial_Node_Data(this)
            spatial_nodes = linspace(0,1,this.n_y)';
        end

        function [time_nodes] = Load_Time_Node_Data(this)
            time_nodes = linspace(0,1,this.n_t)';
        end

        function this = MD_u_Hyperparameters_transient_multi_state_synthetic_test(data_interface, n_y, n_t, component_id)
            this@MD_u_Hyperparameters(data_interface, true, false, true, component_id);
            this.n_y = n_y;
            this.n_t = n_t;
        end

    end

end
