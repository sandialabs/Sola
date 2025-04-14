classdef MD_Data_Interface_synthetic_test_control_vec < MD_Data_Interface

    properties
        m
        x
        epsilon
    end

    methods (Access = public)

        function [u_opt] = Load_Optimal_u(this)
            u_opt = 1 + this.x;
        end

        function [z_opt] = Load_Optimal_z(this)
            z_opt = [1 ; 1];
        end

        function [Z] = Load_Z_Data(this)
            Z = zeros(2, 2);
            Z(:, 1) = [1 ; 1];
            Z(:, 2) = [2 ; 1];
        end

        function [D] = Load_d_Data(this)
            D = zeros(this.m,2);
            D(:,1) = 2*this.epsilon + this.epsilon * this.x;
            D(:,2) = 2*this.epsilon + this.epsilon * this.x;
        end

    end

    methods

        function this = MD_Data_Interface_synthetic_test_control_vec(m)
            this.m = m;
            this.x = linspace(0, 1, m)';
            this.epsilon = 0.1;
        end

    end

end
