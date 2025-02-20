classdef MD_Hyperparameters_Transient_ADR_2D < MD_Hyperparameters

    properties
        x
        y
        t
    end

    methods (Access = public)

        function [spatial_nodes] = Load_Spatial_Node_Data(this)
            spatial_nodes = [this.x, this.y];
        end

        function [time_nodes] = Load_Time_Node_Data(this)
            time_nodes = this.t;
        end

        function this = MD_Hyperparameters_Transient_ADR_2D(data_interface, solver, t)
            this@MD_Hyperparameters(data_interface, true);
            this.x = solver.x;
            this.y = solver.y;
            this.t = t;
        end

    end

end
