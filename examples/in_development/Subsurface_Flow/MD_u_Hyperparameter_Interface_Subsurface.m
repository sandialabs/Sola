%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef MD_u_Hyperparameter_Interface_Subsurface < MD_u_Hyperparameter_Interface
    % Hyperparameter data interface for the subsurface pressure-discrepancy
    % prior.  The state discrepancy is a stationary one-dimensional spatial
    % field on the same nodes as the pressure variable.

    properties
        x
    end

    methods (Access = public)

        function [nodes] = Load_Spatial_Node_Data(this)
            nodes = cell(1, 1);
            nodes{1} = this.x;
        end

        function this = MD_u_Hyperparameter_Interface_Subsurface(x, center_data)
            if nargin < 2
                center_data = false;
            end
            this@MD_u_Hyperparameter_Interface(false, center_data);
            this.x = x;
        end

    end

end
