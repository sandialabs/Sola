classdef Adv_Diff_Constraint < Dynamic_Constraint

    properties
        A
        M
        x
    end

    methods (Access = public)

        function [f, f_y, f_z] = f(this, y, z, t)
            % Evaluate the RHS and its derivatives.
            f = linsolve(this.M, -this.A * y + (1.e2) * this.M * z);
            f_y = linsolve(this.M, -this.A);
            f_z = (1.e2) * eye(this.n_y);
        end

        function [h, h_z] = h(this, z)
            h = zeros(this.n_y, 1);
            h_z = zeros(this.n_y, this.n_z);
        end

        function [Mv] = f_yy_Apply(this, v, y, z, t, lambda)
            num_vecs = size(v, 2);
            Mv = zeros(this.n_y, num_vecs);
        end

        function [Mv] = f_yz_Apply(this, v, y, z, t, lambda)
            num_vecs = size(v, 2);
            Mv = zeros(this.n_y, num_vecs);
        end

        function [Mv] = f_zy_Apply(this, v, y, z, t, lambda)
            num_vecs = size(v, 2);
            Mv = zeros(this.n_z, num_vecs);
        end

        function [Mv] = f_zz_Apply(this, v, y, z, t, lambda)
            num_vecs = size(v, 2);
            Mv = zeros(this.n_z, num_vecs);
        end

        function [Mv] = h_zz_Apply(this, v, z, lambda)
            num_vecs = size(v, 2);
            Mv = zeros(this.n_z, num_vecs);
        end

    end

    methods (Access = public)

        function this = Adv_Diff_Constraint(m, n, T, N)
            this = this@Dynamic_Constraint(m, n, T, N);

            Pe = 5; % Peclet number

            % Spatial domain
            this.x = linspace(0, 1, m)';
            dx = this.x(2) - this.x(1);

            % Mass matrix
            M = diag(4 * ones(1, m)) + diag(ones(1, m - 1), 1) + diag(ones(1, m - 1), -1);
            M(1, 1) = .5 * M(1, 1);
            M(end, end) = .5 * M(end, end);
            M = (1 / 6) * dx * M;

            % Stiffness matrix (diffusion)
            S = diag(2 * ones(1, m)) + (-1) * diag(ones(1, m - 1), 1) + (-1) * diag(ones(1, m - 1), -1);
            S(1, 1) = .5 * S(1, 1);
            S(end, end) = .5 * S(end, end);
            S = (1 / dx) * S;

            % Viscosity matrix (advection)
            V = diag(0 * ones(1, m)) + (1 / 2) * diag(ones(1, m - 1), 1) + (-1 / 2) * diag(ones(1, m - 1), -1);
            V(1, 1) = -1 / 2;
            V(end, end) = 1 / 2;

            % Discretized PDE: M f(y,z) = -Ay + Mz, A = S + Pe V.
            A = S + Pe * V;

            % Store properties.
            this.A = A;
            this.M = M;
        end

    end

end
