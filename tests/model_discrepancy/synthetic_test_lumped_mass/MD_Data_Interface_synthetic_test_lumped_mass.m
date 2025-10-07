classdef MD_Data_Interface_synthetic_test_lumped_mass < MD_Data_Interface

    properties
        m
        x
    end

    methods (Access = public)

        function [u_opt] = Load_Optimal_u(this)
            u_opt = (1 + this.x).^3;
        end

        function [z_opt] = Load_Optimal_z(this)
            z_opt = 1 + this.x;
        end

        function [Z] = Load_Z_Data(this)
            Z = zeros(this.m, 2);
            Z(:, 1) = 1 + this.x;
            Z(:, 2) = this.x + this.x.^2;
        end

        function [D] = Load_d_Data(this)
            Z = this.Load_Z_Data();
            D = .2 * (Z.^3);
        end

    end

    methods

        function this = MD_Data_Interface_synthetic_test_lumped_mass(m)
            this.m = m;
            this.x = linspace(0, 1, m)';
        end

    end

end
