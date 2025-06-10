classdef MD_Data_Interface_synthetic_test_transient < MD_Data_Interface

    properties
        n_y
        n_t
        x
        t
        c_low
        c_high
    end

    methods (Access = public)

        function [u_opt] = Load_Optimal_u(this)
            u_opt = zeros(this.n_y * this.n_t, 1);
            u_opt(1:this.n_y) = (1 + this.x).^3;
            for k = 1:(this.n_t - 1)
                u_opt((k * this.n_y + 1):((k + 1) * this.n_y)) = this.c_low * u_opt(((k - 1) * this.n_y + 1):(k * this.n_y));
            end
        end

        function [z_opt] = Load_Optimal_z(this)
            z_opt = 1 + this.x;
        end

        function [Z] = Load_Z_Data(this)
            Z = zeros(this.n_y, 2);
            Z(:, 1) = 1 + this.x;
            Z(:, 2) = this.x + this.x.^2;
        end

        function [D] = Load_d_Data(this)
            D = zeros(this.n_y * this.n_t, 2);
            for k = 1:this.n_t
                D(((k - 1) * this.n_y + 1):(k * this.n_y), :) = (this.c_high^(k - 1) - this.c_low^(k - 1)) * this.Load_Z_Data().^3;
            end
        end

    end

    methods

        function this = MD_Data_Interface_synthetic_test_transient(n_y, n_t, T, c_low, c_high)
            this.n_y = n_y;
            this.n_t = n_t;
            this.x = linspace(0, 1, n_y)';
            this.t = linspace(0, T, n_t)';
            this.c_low = c_low;
            this.c_high = c_high;
        end

    end

end
