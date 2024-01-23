classdef Thermal_Constraint_AD < Dynamic_Constraint_AD

    properties
        x
        M
        Minv
        S
        forcing
        nodes_to_coll_points
    end

    methods (Access = public)

        function [f] = f_AD(this, y, z, t)
            D = this.Assembly(z);
            f = -this.Minv * D * y + (10^3) * this.forcing(this.x, t);
        end

        function [h] = h_AD(this, z)
            h = ones(this.n_y, 1);
        end

    end

    methods (Access = public)

        function [val] = Forcing_Function(this, x, t)
            val = this.forcing(x, t);
        end

        function [D] = Assembly(this, z)
            h = this.x(2) - this.x(1);
            z_coll = this.nodes_to_coll_points * z;
            diff_x = reshape(z_coll, 2, this.n_y - 1);
            s = sum(diff_x, 1);
            D = diag(([0, s] + [s, 0]) * (1 / h) / 2) + (-1) * diag(s, 1) * (1 / h) / 2 + (-1) * diag(s, -1) * (1 / h) / 2;
        end

        function this = Thermal_Constraint_AD(m, n, T, N)
            this = this@Dynamic_Constraint_AD(m, n, T, N);

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

            this.Minv = inv(this.M);

            S = diag(2 * ones(1, m)) + (-1) * diag(ones(1, m - 1), 1) + (-1) * diag(ones(1, m - 1), -1);
            S(1, 1) = .5 * S(1, 1);
            S(end, end) = .5 * S(end, end);
            S = (1 / h) * S;
            this.S = S;

            this.forcing = @(x, t) this.x.^2;

        end

    end
end
