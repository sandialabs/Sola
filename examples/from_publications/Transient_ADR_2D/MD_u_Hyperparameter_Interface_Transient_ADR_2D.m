classdef MD_u_Hyperparameter_Interface_Transient_ADR_2D < MD_u_Hyperparameter_Interface

    properties
        x
        y
        t
    end

    methods (Access = public)

        function [spatial_nodes] = Load_Spatial_Node_Data(this)
            spatial_nodes = cell(2, 1);
            spatial_nodes{1} = [this.x, this.y];
            spatial_nodes{2} = [this.x, this.y];
        end

        function [time_nodes] = Load_Time_Node_Data(this)
            time_nodes = this.t;
        end

        function this = MD_u_Hyperparameter_Interface_Transient_ADR_2D(is_transient, center_data, adapt_time_variance, component_id, solver, t)
            this@MD_u_Hyperparameter_Interface(is_transient, center_data, adapt_time_variance, component_id);
            this.x = solver.x;
            this.y = solver.y;
            this.t = t;
        end

    end

end
