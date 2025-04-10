classdef MD_u_Hyperparameter_Interface_Py < MD_u_Hyperparameter_Interface

    properties
        u_hyperparam_interface_py
    end

    methods

        function [spatial_nodes] = Load_Spatial_Node_Data(this)
            spatial_nodes = this.u_hyperparam_interface_py.Load_Spatial_Node_Data();
            spatial_nodes = cell(spatial_nodes);
            for k = 1:length(spatial_nodes)
                spatial_nodes{k} = double(spatial_nodes{k});
                if size(spatial_nodes{k}, 2) > 3
                    spatial_nodes{k} = spatial_nodes{k}';
                end
            end
        end

        function [time_nodes] = Load_Time_Node_Data(this)
            time_nodes = this.u_hyperparam_interface_py.Load_Time_Node_Data();
            time_nodes = double(time_nodes);
            if size(time_nodes, 2) > 1
                time_nodes = time_nodes';
            end
        end

        function this = MD_u_Hyperparameter_Interface_Py(u_hyperparam_interface_py, is_transient, center_data, adapt_time_variance, component_id)
            arguments
                u_hyperparam_interface_py 
                is_transient {boolean}
                center_data = false
                adapt_time_variance = false
                component_id = 1
            end
            this@MD_u_Hyperparameter_Interface(is_transient, center_data, adapt_time_variance, component_id);
            this.u_hyperparam_interface_py = u_hyperparam_interface_py;
        end
    end
end
