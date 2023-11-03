classdef Thermal_Constraint < Constraint

    properties
        m
        x
        dirichlet_bc
        M
        S
        forcing
        B
        nodes_to_coll_points
    end

    methods (Access = public)

        function [u] = State_Solve(this, z)
            D = this.Assembly(z);
            f = this.B * this.forcing + this.dirichlet_bc;
            u = linsolve(D, f);
        end

        function [c] = c(this, u, z)
            D = this.Assembly(z);
            f = this.B * this.forcing + this.dirichlet_bc;
            c = D * u - f;
        end

        function [Mv] = c_u_Transpose_Inverse_Apply(this, v, u, z)
            D = this.Assembly(z);
            Mv = linsolve(D', v);
        end

        function [Mv] = c_z_Transpose_Apply(this, v, u, z)
            D_diff = this.Assembly_z_Jacobian(u);
            Mv = D_diff' * v;
        end

        function [Mv] = c_u_Inverse_Apply(this, v, u, z)
            D = this.Assembly(z);
            Mv = linsolve(D, v);
        end

        function [Mv] = c_z_Apply(this, v, u, z)
            D_diff = this.Assembly_z_Jacobian(u);
            Mv = D_diff * v;
        end

        function [Mv] = c_uu_Apply(this, v, u, z, lambda)
            Mv = zeros(this.m, 1);
        end

        function [Mv] = c_uz_Apply(this, v, u, z, lambda)
            Mv = zeros(this.m, size(v, 2));
            for k = 1:size(v, 2)
                D = this.Assembly(z);
                D_pert = this.Assembly(z + v(:, k));
                Mv(:, k) = D_pert' * lambda - D' * lambda;
            end
        end

        function [Mv] = c_zu_Apply(this, v, u, z, lambda)
            Mv = zeros(this.m, size(v, 2));
            for k = 1:size(v, 2)
                D_diff = this.Assembly_z_Jacobian(u);
                D_diff_pert = this.Assembly_z_Jacobian(u + v(:, k));
                Mv(:, k) = D_diff_pert' * lambda - D_diff' * lambda;
            end
        end

        function [Mv] = c_zz_Apply(this, v, u, z, lambda)
            Mv = zeros(this.m, 1);
        end

    end

    methods (Access = public)

        function [D] = Assembly(this, z)
            h = this.x(2) - this.x(1);
            z_coll = this.nodes_to_coll_points * z;
            diff_x = reshape(z_coll, 2, this.m - 1);
            s = sum(diff_x, 1);
            D = diag(([0, s] + [s, 0]) * (1 / h) / 2) + (-1) * diag(s, 1) * (1 / h) / 2 + (-1) * diag(s, -1) * (1 / h) / 2;
            D(1, :) = 0 * D(1, :);
            D(this.m, :) = 0 * D(this.m, :);
            D(1, 1) = 1;
            D(this.m, this.m) = 1;
        end

        function [D_diff] = Assembly_z_Jacobian(this, u)
            h = this.x(2) - this.x(1);
            up = (u(2:end) - u(1:end - 1)) / h;
            d = [-(1 / 2) * up(1); (1 / 2) * up(1:end - 1) - (1 / 2) * up(2:end); (1 / 2) * up(end)];
            D_diff = diag(d) + (-1 / 2) * diag(up, 1) + (1 / 2) * diag(up, -1);
            D_diff(1, :) = 0 * D_diff(1, :);
            D_diff(this.m, :) = 0 * D_diff(this.m, :);
        end

        function this = Thermal_Constraint(m)
            this = this@Constraint();

            this.m = m;
            this.x = linspace(0, 1, m)';

            h = this.x(2) - this.x(1);

            coll_points = zeros(2 * (m - 1), 1);
            for k = 1:(m - 1)
                map_to_coll = (1:2)' + 2 * (k - 1);
                coll_points(map_to_coll) = this.x(k) + h * ((1 / sqrt(3)) * [-1; 1] + 1) / 2;
            end

            nodes_to_coll_points = zeros(2 * (m - 1), m);
            for k = 1:(m - 1)
                map_to_coll = (1:2)' + 2 * (k - 1);
                nodes_to_coll_points(map_to_coll(1), k) = (coll_points(map_to_coll(1)) - this.x(k + 1)) / (this.x(k) - this.x(k + 1));
                nodes_to_coll_points(map_to_coll(1), k + 1) = (coll_points(map_to_coll(1)) - this.x(k)) / (this.x(k + 1) - this.x(k));
                nodes_to_coll_points(map_to_coll(2), k) = (coll_points(map_to_coll(2)) - this.x(k + 1)) / (this.x(k) - this.x(k + 1));
                nodes_to_coll_points(map_to_coll(2), k + 1) = (coll_points(map_to_coll(2)) - this.x(k)) / (this.x(k + 1) - this.x(k));
            end
            this.nodes_to_coll_points = nodes_to_coll_points;

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

            this.dirichlet_bc = zeros(m, 1);
            this.dirichlet_bc(1) = 0;
            this.dirichlet_bc(m) = 0;

            this.forcing = this.x.^2;

            B = (10^3) * this.M;
            B(1, :) = 0 * B(1, :);
            B(end, :) = 0 * B(end, :);
            this.B = B;

        end

        function [diffs] = Jacobian_Test(this)
            z = randn(this.m, 1);
            u = rand(this.m, 1);
            h = 10.^(-2:-1:-6);
            p = length(h);

            D_diff = this.Assembly_z_Jacobian(u);
            v = randn(this.m, 1);
            D_diff_v = D_diff * v;
            D = this.Assembly(z);

            diffs = zeros(p, 1);
            for k = 1:p
                D_pert = this.Assembly(z + h(k) * v);
                FD = (D_pert * u - D * u) / h(k);
                diffs(k) = norm(D_diff_v - FD);
            end
            disp('Jacobian finite difference check');
            for k = 1:p
                disp(['h = ', num2str(h(k)), ' and error = ', num2str(diffs(k))]);
            end
            disp(' ');
        end

    end
end
