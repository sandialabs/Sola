%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef MD_z_Hyperparameter_Interface_Subsurface < MD_z_Hyperparameter_Interface
    % Hyperparameter data interface for the subsurface control prior.  The
    % injection/production control is represented as a one-dimensional
    % stationary spatial field on the same mesh as the pressure state.

    properties
        x
        con
    end

    methods (Access = public)

        function [nodes] = Load_Spatial_Node_Data(this)
            nodes = this.x;
        end

        function [u] = State_Solve(this, z)
            u = this.con.State_Solve(z);
        end

        function this = MD_z_Hyperparameter_Interface_Subsurface(x, con, num_state_solves)
            if nargin < 3
                num_state_solves = 0;
            end
            this@MD_z_Hyperparameter_Interface('spatial field', num_state_solves);
            this.x = x;
            this.con = con;
        end

    end

end
