classdef MD_Opt_Prob_Interface_synthetic_test_transient < MD_Opt_Prob_Interface

    properties
        n_y
        n_t
        c_low
        x
        t
        M
    end

    methods (Access = public)

        function [z_out] = Apply_Solution_Operator_z_Jacobian_Transpose(this, u_in, z)
            J = zeros(this.n_y * this.n_t, this.n_y);
            for k = 1:this.n_t
                I = (1:this.n_y) + (k - 1) * this.n_y;
                J(I, :) = (this.c_low^(k - 1)) * 3 * diag(z.^2);
            end
            z_out = J' * u_in;
        end

        function [z_out] = Apply_RS_Hessian(this, z_in, z)
            tmp1 = (this.c_low^(this.n_t - 1)) * 3 * diag(z.^2) * z_in;
            tmp2 = this.M * tmp1;
            z_out = (this.c_low^(this.n_t - 1)) * 3 * diag(z.^2) * tmp2;
        end

        function [grad_u] = Misfit_Gradient(this, u, z)
            grad_u = 0 * u;
            grad_u(((this.n_t - 1) * this.n_y + 1):end) = this.M * (u(((this.n_t - 1) * this.n_y + 1):end) - (this.c_low^(this.n_t - 1)) * (1 + this.x).^3);
        end

        function [u_out] = Apply_Misfit_Hessian(this, u_in, u, z)
            u_out = 0 * u_in;
            u_out(((this.n_t - 1) * this.n_y + 1):end, :) = this.M * u_in(((this.n_t - 1) * this.n_y + 1):end, :);
        end

    end

    methods

        function this = MD_Opt_Prob_Interface_synthetic_test_transient(n_y, n_t, T, c_low)
            this.n_y = n_y;
            this.n_t = n_t;
            this.c_low = c_low;
            this.x = linspace(0, 1, n_y)';
            this.t = linspace(0, T, n_t)';

            h = this.x(2) - this.x(1);
            M = diag(4 * ones(1, n_y)) + diag(ones(1, n_y - 1), 1) + diag(ones(1, n_y - 1), -1);
            M(1, 1) = .5 * M(1, 1);
            M(end, end) = .5 * M(end, end);
            M = (1 / 6) * h * M;
            this.M = M;

        end

    end

end
