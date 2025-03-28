classdef MD_u_Hyperparameter_Interface_Transient_Test_Problem < MD_u_Hyperparameter_Interface

    properties
        x
        t
    end

    methods (Access = public)

        function [nodes] = Load_Spatial_Node_Data(this)
            nodes = cell(1, 1);
            nodes{1} = this.x;
        end

        function [time_nodes] = Load_Time_Node_Data(this)
            time_nodes = this.t;
        end

        function this = MD_u_Hyperparameter_Interface_Transient_Test_Problem(n_y, n_t)
            this@MD_u_Hyperparameter_Interface(true, false, true);
            this.x = linspace(0, 1, n_y)';
            this.t = linspace(0, 1, n_t)';
        end

    end

end
