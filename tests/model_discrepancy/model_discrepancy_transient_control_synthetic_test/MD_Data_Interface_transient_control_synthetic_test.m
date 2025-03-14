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
            tmp = [this.t,2*this.t]';
            z_opt = tmp(:);
        end

        function [Z] = Load_Z_Data(this)
            Z = this.Load_Optimal_z();
        end

        function [D] = Load_d_Data(this)
            D = ones(this.n_y*this.n_t,1);
        end

    end

    methods

        function this = MD_Data_Interface_transient_control_synthetic_test(n_y,n_t)
            this.n_y = n_y;
            this.x = linspace(0, 1, n_y)';
            this.n_t = n_t;
            this.t = linspace(0,1,n_t)';

            J2 = (0:(n_y-1))/(n_y-1);
            tmp = 2*J2'*this.t' + ones(n_y,1)*this.t';
            this.T = tmp(:);
        end

    end

end
