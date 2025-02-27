classdef MD_z_Hyperparameters_Transient_Test_Problem < MD_z_Hyperparameters

    properties
        x
        t
    end

    methods (Access = public)

        function [nodes] = Load_Spatial_Node_Data(this)
            nodes = this.x;
        end

        function this = MD_z_Hyperparameters_Transient_Test_Problem(data_interface, u_prior_interface, n_y)
            this@MD_z_Hyperparameters(data_interface, u_prior_interface);
            this.x = linspace(0, 1, n_y)';
            n_t = length(data_interface.u_opt) / n_y;
            this.t = linspace(0, 1, n_t)';
        end

    end

end
