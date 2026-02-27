classdef MD_Data_Interface_lumped_mass_unit_test < MD_Data_Interface

    properties
        n_y
        n_t
        x
        t
    end

    methods (Access = public)

        function [u_opt] = Load_Optimal_u(this)
            u_opt = (1 + this.x).^3 .* this.t';
            u_opt = u_opt(:);
        end

        function [z_opt] = Load_Optimal_z(this)
            z_opt = 1 + this.x;
        end

        function [Z] = Load_Z_Data(this)
            Z = 1 + this.x;
        end

        function [D] = Load_d_Data(this)
            Z = this.Load_Z_Data();
            D = .2 * Z .* this.t';
            D = D(:);
        end

    end

    methods

        function this = MD_Data_Interface_lumped_mass_unit_test(n_y, n_t)
            this.n_y = n_y;
            this.n_t = n_t;
            this.x = linspace(0, 1, n_y)';
            this.t = linspace(0, 1, n_t)';
        end

    end

end
