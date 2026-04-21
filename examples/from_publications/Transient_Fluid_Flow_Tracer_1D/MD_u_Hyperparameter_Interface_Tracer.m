%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef MD_u_Hyperparameter_Interface_Tracer < MD_u_Hyperparameter_Interface

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

        function this = MD_u_Hyperparameter_Interface_Tracer(n_t, n_y, adapt_time_variance)
            this@MD_u_Hyperparameter_Interface(true, false, adapt_time_variance);
            this.x = linspace(0, 1, n_y)';
            this.t = linspace(0, 1, n_t)';
        end

    end

end
