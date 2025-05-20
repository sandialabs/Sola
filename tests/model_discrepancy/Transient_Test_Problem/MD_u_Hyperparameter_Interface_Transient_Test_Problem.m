classdef MD_u_Hyperparameter_Interface_Transient_Test_Problem < MD_u_Hyperparameter_Interface

    properties
        x
        t
    end

    methods (Access = public)

        function [spatial_nodes] = Load_Spatial_Node_Data(this)
            spatial_nodes = cell(1, 1);
            spatial_nodes{1} = this.x;
        end

        function [time_nodes] = Load_Time_Node_Data(this)
            time_nodes = this.t;
        end

        function this = MD_u_Hyperparameter_Interface_Transient_Test_Problem(x, t, adapt_time_variance)
            this@MD_u_Hyperparameter_Interface(true,false,adapt_time_variance);
            this.x = x;
            this.t = t;
        end

    end

end
