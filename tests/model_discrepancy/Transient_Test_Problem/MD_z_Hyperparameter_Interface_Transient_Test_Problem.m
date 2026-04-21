%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef MD_z_Hyperparameter_Interface_Transient_Test_Problem < MD_z_Hyperparameter_Interface

    properties
        x
        con
    end

    methods (Access = public)

        function [nodes] = Load_Spatial_Node_Data(this)
            nodes = this.x;
        end

        function [u] = State_Solve(this, z)

            u1 = this.con.State_Solve(z(:, 1));
            m = length(u1);
            n = size(z, 2);
            u = zeros(m, n);
            u(:, 1) = u1;
            for k = 2:n
                u(:, k) = this.con.State_Solve(z(:, k));
            end
        end

        function this = MD_z_Hyperparameter_Interface_Transient_Test_Problem(x, con)
            this@MD_z_Hyperparameter_Interface('spatial field', 100);
            this.x = x;
            this.con = con;
        end

    end

end
