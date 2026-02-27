classdef MD_Data_Interface_transient_control_synthetic_test < MD_Data_Interface

    properties
        n_y
        n_t
        x
        t
        T
    end

    methods (Access = public)

        function [u_opt] = Load_Optimal_u(this)
            u_opt = this.T;
        end

        function [z_opt] = Load_Optimal_z(this)
            tmp = [this.t, 2 * this.t]';
            z_opt = tmp(:);
        end

        function [Z] = Load_Z_Data(this)
            Z = this.Load_Optimal_z();
        end

        function [D] = Load_d_Data(this)
            D = ones(this.n_y * this.n_t, 1);
        end

    end

    methods

        function this = MD_Data_Interface_transient_control_synthetic_test(n_y, n_t)
            this.n_y = n_y;
            this.x = linspace(0, 1, n_y)';
            this.n_t = n_t;
            this.t = linspace(0, 1, n_t)';

            J1 = 1 - (0:(n_y - 1)) / (n_y - 1);
            J2 = (0:(n_y - 1)) / (n_y - 1);
            J = kron(eye(n_t), [J1; J2])';

            % Forward model: u(t) = z_1(t) * J1 + z_2(t) * J2

            tmp = [this.t, 2 * this.t]';
            z_opt = tmp(:);
            this.T = J * z_opt;
        end

    end

end
