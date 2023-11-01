classdef Poisson_Constraint < Constraint

    properties
        m
        diff_coeff
        x
        M
        S
        A
        B
    end

    methods (Access = public)

        function [u] = State_Solve(this, z)
            b = this.B * z;
            u = linsolve(this.A, b);
        end

        function [Mv] = c_u_Transpose_Inverse_Apply(this, v, u, z)
            Mv = linsolve(this.A', v);
        end

        function [Mv] = c_z_Transpose_Apply(this, v, u, z)
            Mv = -this.B' * v;
        end

        function [Mv] = c_u_Inverse_Apply(this, v, u, z)
            Mv = linsolve(this.A, v);
        end

        function [Mv] = c_z_Apply(this, v, u, z)
            Mv = -this.B * v;
        end

        function [Mv] = c_uu_Apply(this, v, u, z, lambda)
            Mv = zeros(this.m, 1);
        end

        function [Mv] = c_uz_Apply(this, v, u, z, lambda)
            Mv = zeros(this.m, 1);
        end

        function [Mv] = c_zu_Apply(this, v, u, z, lambda)
            Mv = zeros(this.m, 1);
        end

        function [Mv] = c_zz_Apply(this, v, u, z, lambda)
            Mv = zeros(this.m, 1);
        end

    end

    methods (Access = public)

        function this = Poisson_Constraint(m, diff_coeff)
            this = this@Constraint();
            this.m = m;
            this.diff_coeff = diff_coeff;
            this.x = linspace(0, 1, m)';

            h = this.x(2) - this.x(1);

            M = diag(4 * ones(1, m)) + diag(ones(1, m - 1), 1) + diag(ones(1, m - 1), -1);
            M(1, 1) = .5 * M(1, 1);
            M(end, end) = .5 * M(end, end);
            M = (1 / 6) * h * M;
            this.M = M;

            S = diag(2 * ones(1, m)) + (-1) * diag(ones(1, m - 1), 1) + (-1) * diag(ones(1, m - 1), -1);
            S(1, 1) = .5 * S(1, 1);
            S(end, end) = .5 * S(end, end);
            S = (1 / h) * S;
            this.S = S;

            A = diff_coeff * S;
            A(1, :) = 0 * A(1, :);
            A(end, :) = 0 * A(end, :);
            A(1, 1) = 1;
            A(end, end) = 1;
            this.A = A;

            B = M;
            B(1, :) = 0 * B(1, :);
            B(end, :) = 0 * B(end, :);
            this.B = B;
        end

    end
end
