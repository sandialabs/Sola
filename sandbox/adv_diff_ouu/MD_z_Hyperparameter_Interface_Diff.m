classdef MD_z_Hyperparameter_Interface_Diff < MD_z_Hyperparameter_Interface

    properties
        x
        cons_lofi
    end

    methods (Access = public)

        function [nodes] = Load_Spatial_Node_Data(this)
            nodes = this.x;
        end

        function [u] = State_Solve(this, z)
            N = size(z, 2);
            m = length(this.x);
            n_r = length(this.cons_lofi);
            tmp = zeros(m, N, n_r);
            for k = 1:n_r
                tmp(:, :, k) = this.cons_lofi{k}.State_Solve(z);
            end
            u = mean(tmp, 3);
        end

        function this = MD_z_Hyperparameter_Interface_Diff(num_state_solves, x, cons_lofi)
            this@MD_z_Hyperparameter_Interface('spatial field', num_state_solves);
            this.x = x;
            this.cons_lofi = cons_lofi;
        end

    end

end
