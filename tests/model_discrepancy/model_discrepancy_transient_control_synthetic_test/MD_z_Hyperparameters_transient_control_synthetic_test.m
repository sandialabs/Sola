classdef MD_z_Hyperparameters_transient_control_synthetic_test < MD_z_Hyperparameters

    properties
        t
        opt_prob_interface
    end

    methods (Access = public)

        function [time_nodes] = Load_Time_Node_Data(this)
            time_nodes = this.t;
        end

        function [u] = State_Solve(this,z)
            u = this.opt_prob_interface.J*z;
        end

        function this = MD_z_Hyperparameters_transient_control_synthetic_test(data_interface, u_prior_interface, num_state_solves, t, opt_prob_interface)
            this@MD_z_Hyperparameters(data_interface, u_prior_interface,'transient vector',num_state_solves);
            this.t = t;
            this.opt_prob_interface = opt_prob_interface;
        end

    end

end
