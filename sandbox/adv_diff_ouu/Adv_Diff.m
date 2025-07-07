classdef Adv_Diff < Parametric_Constraint

    properties
        m
        diff_coeff
        vel_coeff
        x
        M
        S
        V
        A
        B
        c
    end

    methods (Access = public)

        function [u] = Parametric_State_Solve(this, z, theta)
            if norm(theta - this.theta_current) ~= 0
                this.Assemble(theta);
            end
            u = linsolve(this.A, this.B * z + this.c);
        end

        function [Mv] = Parametric_c_u_Transpose_Inverse_Apply(this, v, u, z, theta)
            if norm(theta - this.theta_current) ~= 0
                this.Assemble(theta);
            end
            Mv = linsolve(this.A', v);
        end

        function [Mv] = Parametric_c_z_Transpose_Apply(this, v, u, z, theta)
            Mv = -this.B' * v;
        end

        function [Mv] = Parametric_c_u_Inverse_Apply(this, v, u, z, theta)
            if norm(theta - this.theta_current) ~= 0
                this.Assemble(theta);
            end
            Mv = linsolve(this.A, v);
        end

        function [Mv] = Parametric_c_z_Apply(this, v, u, z, theta)
            Mv = -this.B * v;
        end

        function [Mv] = Parametric_c_uu_Apply(this, v, u, z, lambda, theta)
            Mv = zeros(this.m, 1);
        end

        function [Mv] = Parametric_c_uz_Apply(this, v, u, z, lambda, theta)
            Mv = zeros(this.m, 1);
        end

        function [Mv] = Parametric_c_zu_Apply(this, v, u, z, lambda, theta)
            Mv = zeros(this.m, 1);
        end

        function [Mv] = Parametric_c_zz_Apply(this, v, u, z, lambda, theta)
            Mv = zeros(this.m, 1);
        end

    end

    methods (Access = public)

        function [] = Assemble(this, theta)
            this.diff_coeff = theta;
            A = this.diff_coeff * this.S + this.vel_coeff * this.V;
            A(1, :) = 0;
            A(1, 1) = 1;
            this.A = A;
        end

        function this = Adv_Diff(m, vel_coeff, theta)
            this = this@Parametric_Constraint(theta);
            this.m = m;
            this.vel_coeff = vel_coeff;
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

            V = diag(0 * ones(1, m)) + (1 / 2) * diag(ones(1, m - 1), 1) + (-1 / 2) * diag(ones(1, m - 1), -1);
            V(1, 1) = -1 / 2;
            V(end, end) = 1 / 2;
            this.V = V;

            B = (10^2) * this.M;
            B(1, :) = 0;
            this.B = B;

            c = zeros(this.m, 1);
            c(1) = 1;
            this.c = c;

            this.Assemble(theta);
        end

    end
end
