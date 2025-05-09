classdef MD_Data_Interface_transient_multi_state_synthetic < MD_Data_Interface

    properties
        n_y
        n_t
        x
        c_low
        c_high
    end

    methods (Access = public)

        function [u_opt] = Load_Optimal_u(this)
            u1 = kron(this.c_low.^(0:(this.n_t - 1))', (1 + this.x).^3);
            u2 = this.c_low * u1;
            u_tmp = [reshape(u1, 50, 10); reshape(u2, 50, 10)];
            u_opt = u_tmp(:);
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
            Z = this.Load_Z_Data();

            D = zeros(2 * this.n_y * this.n_t, 2);
            for k = 1:2
                u1 = kron(this.c_low.^(0:(this.n_t - 1))', (Z(:,k)).^3);
                u2 = this.c_low * u1;
                u_tmp = [reshape(u1, 50, 10); reshape(u2, 50, 10)];
                u_lofi = u_tmp(:);

                u1 = kron(this.c_high.^(0:(this.n_t - 1))', (Z(:,k)).^3);
                u2 = this.c_high * u1;
                u_tmp = [reshape(u1, 50, 10); reshape(u2, 50, 10)];
                u_hifi = u_tmp(:);

                D(:, k) = u_hifi - u_lofi;
            end

        end

        % Defaults to return all state elements
        % Overload function to extract component i from the state
        % Returns a vector of integers I index elements of component i
        function [I] = Separate_State_Components(this, i)
            if i == 1
                I = zeros(this.n_y * this.n_t, 1);
                for k = 1:this.n_t
                    I(((k - 1) * this.n_y + 1):(k * this.n_y)) = (1:this.n_y) + 2 * (k - 1) * this.n_y;
                end
            else
                I = zeros(this.n_y * this.n_t, 1);
                for k = 1:this.n_t
                    I(((k - 1) * this.n_y + 1):(k * this.n_y)) = ((this.n_y + 1):(2 * this.n_y)) + 2 * (k - 1) * this.n_y;
                end
            end
        end

    end

    methods

        function this = MD_Data_Interface_transient_multi_state_synthetic (n_y, n_t, c_low, c_high)
            this.n_y = n_y;
            this.n_t = n_t;
            this.x = linspace(0, 1, n_y)';
            this.c_low = c_low;
            this.c_high = c_high;
        end

    end

end
