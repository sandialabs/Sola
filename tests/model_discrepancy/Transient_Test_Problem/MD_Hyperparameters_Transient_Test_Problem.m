classdef MD_Hyperparameters_Transient_Test_Problem < MD_Hyperparameters

    properties
        x
        t
    end

    methods (Access = public)

        function [spatial_nodes] = Load_Spatial_Node_Data(this)
            spatial_nodes = this.x;
        end

        function [time_nodes] = Load_Time_Node_Data(this)
            time_nodes = this.t;
        end

        function this = MD_Hyperparameters_Transient_Test_Problem(data_interface, x, t)
            this@MD_Hyperparameters(data_interface, true);
            this.x = x;
            this.t = t;
        end

    end

end
