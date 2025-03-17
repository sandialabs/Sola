classdef MD_z_Hyperparameter_Interface_Transient_Test_Problem < MD_z_Hyperparameter_Interface

    properties
        x
        t
        con
    end

    methods (Access = public)

        function [nodes] = Load_Spatial_Node_Data(this)
            nodes = this.x;
        end

        function [u] = State_Solve(this,z)
            N = size(z,2);
            u = zeros(length(this.x)*length(this.t),N);
            for k = 1:N
                u(:,k) = this.con.State_Solve(z(:,k));
            end
        end

        function this = MD_z_Hyperparameter_Interface_Transient_Test_Problem(num_state_solves, con, n_y, n_t)
            this@MD_z_Hyperparameter_Interface('spatial field',num_state_solves);
            this.x = linspace(0, 1, n_y)';
            this.t = linspace(0, 1, n_t)';
            this.con = con;
        end

    end

end
