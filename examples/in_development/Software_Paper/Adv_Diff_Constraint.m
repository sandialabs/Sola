classdef Adv_Diff_Constraint < Constraint
    % Advection-diffusion system, discretized with finite elements and
    % integrated with implicit Euler time stepping.

    properties
        m
        diff_coeff  % Diffusion constant :math:`\kappa`
        vel_coeff  % Advection / velocity constant :math:`v`
        robin_coeff  % Robin boundary condition constant :math:`h`
        x  % Spatial domain
        M  % Mass matrix
        S  % Stiffness matrix
        V  % Velocity matrix
        robin_bc  % Boundary condition matrix
        A  % Au = M z.
    end

    methods (Access = public)

        function [u] = State_Solve(this, z)
            u = linsolve(this.A, this.M * z);
        end

        function [Mv] = c_u_Transpose_Inverse_Apply(this, v, u, z)
            Mv = linsolve(this.A', v);
        end

        function [Mv] = c_z_Transpose_Apply(this, v, u, z)
            Mv = -this.M' * v;
        end

        function [Mv] = c_u_Inverse_Apply(this, v, u, z)
            Mv = linsolve(this.A, v);
        end

        function [Mv] = c_z_Apply(this, v, u, z)
            Mv = -this.M * v;
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

        function this = Adv_Diff_Constraint(m, diff_coeff, vel_coeff, robin_coeff)
            this = this@Constraint();
            this.m = m;
            this.diff_coeff = diff_coeff;
            this.vel_coeff = vel_coeff;
            this.robin_coeff = robin_coeff;
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

            robin_bc = zeros(m, m);
            robin_bc(1, 1) = 1;
            robin_bc(end, end) = 1;
            this.robin_bc = robin_bc;

            this.A = this.diff_coeff * this.S + this.vel_coeff * this.V + this.robin_coeff * this.robin_bc;
        end

    end
end
