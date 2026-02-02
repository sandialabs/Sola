classdef MD_z_Hyperparameter_Interface_Tracer < MD_z_Hyperparameter_Interface

    properties
        x
        y
        con_lofi
    end

    methods (Access = public)

        function [nodes] = Load_Spatial_Node_Data(this)
            nodes = this.x;
        end

        % function [u] = State_Solve(this, z)
        %     u = this.con_lofi.State_Solve(z);
        % end

        function this = MD_z_Hyperparameter_Interface_Tracer(x)
            this@MD_z_Hyperparameter_Interface('spatial field');
            this.x = x;
            % this.con_lofi = con_lofi;
        end

    end

end

% function [nodes] = Load_Spatial_Node_Data(this)
%             nodes = this.x;
%         end

%         function this = MD_z_Hyperparameter_Interface_Diff(x)
%             this@MD_z_Hyperparameter_Interface('spatial field');
%             this.x = x;
%         end
