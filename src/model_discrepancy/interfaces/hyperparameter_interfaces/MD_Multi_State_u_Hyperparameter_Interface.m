classdef MD_Multi_State_u_Hyperparameter_Interface < MD_u_Hyperparameter_Interface

    properties
        u_hyperparam_interface_cell
    end

    methods

        function [spatial_nodes] = Load_Spatial_Node_Data(this)
            spatial_nodes = this.u_hyperparam_interface_cell{1}.Load_Spatial_Node_Data();
        end

        function [time_nodes] = Load_Time_Node_Data(this)
            time_nodes = this.u_hyperparam_interface_cell{1}.Load_Time_Node_Data();
        end

        function this = MD_Multi_State_u_Hyperparameter_Interface(u_hyperparam_interface_cell)
            this@MD_u_Hyperparameter_Interface(u_hyperparam_interface_cell{1}.is_transient, u_hyperparam_interface_cell{1}.center_data, u_hyperparam_interface_cell{1}.adapt_time_variance, 0)
            this.u_hyperparam_interface_cell = u_hyperparam_interface_cell;
        end

    end

end
