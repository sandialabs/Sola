classdef MD_u_Hyperparameters_Transient_Test_Problem < MD_u_Hyperparameters

    properties
        x
        t
    end

    methods (Access = public)

        function [nodes] = Load_Spatial_Node_Data(this)
            nodes = this.x;
        end

        function [time_nodes] = Load_Time_Node_Data(this)
            time_nodes = this.t;
        end

        function this = MD_u_Hyperparameters_Transient_Test_Problem(data_interface, n_y)
            this@MD_u_Hyperparameters(data_interface, true, false, true);
            this.x = linspace(0, 1, n_y)';
            n_t = length(data_interface.u_opt) / n_y;
            this.t = linspace(0, 1, n_t)';
        end

    end

end
