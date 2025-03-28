classdef MD_z_Hyperparameters_Transient_ADR_2D < MD_z_Hyperparameters

    properties
        t
    end

    methods (Access = public)

        function [time_nodes] = Load_Time_Node_Data(this)
            time_nodes = this.t;
        end

        function this = MD_z_Hyperparameters_Transient_ADR_2D(data_interface, u_prior_interface, n_t)
            this@MD_z_Hyperparameters(data_interface, u_prior_interface, 'transient vector');
            this.t = linspace(0, 1, n_t)';
            this.t = this.t(1:end - 1);
        end

    end

end
