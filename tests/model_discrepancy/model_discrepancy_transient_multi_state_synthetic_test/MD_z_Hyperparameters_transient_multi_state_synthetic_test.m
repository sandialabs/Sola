classdef MD_z_Hyperparameters_transient_multi_state_synthetic_test < MD_z_Hyperparameters

    properties
        n_y
        n_t
    end

    methods (Access = public)

        function [spatial_nodes] = Load_Spatial_Node_Data(this)
            spatial_nodes = linspace(0,1,this.n_y)';
        end

        function [u] = State_Solve(this,z)
            p = size(z,2);
            u = zeros(2*this.n_y*this.n_t,p);
            for k = 1:p
                u1 = kron(this.data_interface.c_low.^(0:(this.n_t-1))',z(:,k).^3);
                u2 = this.data_interface.c_low * u1;
                u_tmp = [u1 , u2];
                u(:,k) = u_tmp(:);
            end
        end

        function this = MD_z_Hyperparameters_transient_multi_state_synthetic_test(data_interface, u_prior_interface, num_state_solves, n_y, n_t)
            this@MD_z_Hyperparameters(data_interface, u_prior_interface, num_state_solves);
            this.n_y = n_y;
            this.n_t = n_t;
        end

    end

end
