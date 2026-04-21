%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Adv_Diff < Constraint

    properties
        m
        T
        diff_coeff
        vel_coeff
        xi
        x
        M
        S
        V
        A
        B
        c
    end

    methods (Access = public)

        function [u] = State_Solve(this, z)
            u = linsolve(this.A, this.B * z + this.c);
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

        function [Mv] = Objective_uu_Apply(this, v, u, z)
            Mv = this.M * v;
        end

        function [Mv] = Objective_uz_Apply(this, v, u, z)
            Mv = zeros(this.m, 1);
        end

        function [Mv] = Objective_zu_Apply(this, v, u, z)
            Mv = zeros(this.m, 1);
        end

        function [Mv] = Objective_zz_Apply(this, v, u, z)
            Mv = this.reg_coeff * this.reg_mat * v;
        end

    end

    methods (Access = public)

        function this = Adv_Diff(m, vel_coeff, xi)
            this = this@Constraint();
            this.m = m;
            this.xi = xi;
            this.diff_coeff = xi;
            this.vel_coeff = vel_coeff;
            this.x = linspace(0, 1, m)';
            this.T = 50 - 60 * (this.x - 0.5).^2;

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

            A = this.diff_coeff * this.S + this.vel_coeff * this.V;
            A(1, :) = 0;
            A(1, 1) = 1;
            this.A = A;

            B = (10^2) * this.M;
            B(1, :) = 0;
            this.B = B;

            c = zeros(this.m, 1);
            c(1) = 1;
            this.c = c;
        end

    end
end
