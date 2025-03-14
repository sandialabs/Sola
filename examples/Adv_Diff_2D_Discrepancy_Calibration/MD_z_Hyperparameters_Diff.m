classdef MD_z_Hyperparameters_Diff < MD_z_Hyperparameters

    properties
        x
        y
        con_lofi
    end

    methods (Access = public)

        function [nodes] = Load_Spatial_Node_Data(this)
            nodes = [this.x, this.y];
        end

        function [u] = State_Solve(this,z)
            u = this.con_lofi.State_Solve(z);
        end

        function this = MD_z_Hyperparameters_Diff(data_interface, u_prior_interface, num_state_solves, x, y, con_lofi)
            this@MD_z_Hyperparameters(data_interface, u_prior_interface, 'spatial field', num_state_solves);
            this.x = x;
            this.y = y;
            this.con_lofi = con_lofi;
        end

    end

end
