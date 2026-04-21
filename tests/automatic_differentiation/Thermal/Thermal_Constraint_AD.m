%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Thermal_Constraint_AD < Constraint_AD

    properties
        n_y
        x
        dirichlet_bc
        M
        S
        forcing
        B
        nodes_to_coll_points
    end

    methods (Access = public)

        function [c] = c_AD(this, u, z)
            D = this.Assembly(z);
            f = this.B * this.forcing + this.dirichlet_bc;
            c = D * u - f;
        end

    end

    methods (Access = public)

        function [D] = Assembly(this, z)
            h = this.x(2) - this.x(1);
            z_coll = this.nodes_to_coll_points * z;
            diff_x = reshape(z_coll, 2, this.n_y - 1);
            s = sum(diff_x, 1);
            D = diag(([0, s] + [s, 0]) * (1 / h) / 2) + (-1) * diag(s, 1) * (1 / h) / 2 + (-1) * diag(s, -1) * (1 / h) / 2;
            D(1, :) = double(1:this.n_y == 1);
            D(this.n_y, :) = double(1:this.n_y == this.n_y);
        end

        function this = Thermal_Constraint_AD(n_y)
            this = this@Constraint_AD(n_y, n_y);

            this.n_y = n_y;
            this.x = linspace(0, 1, n_y)';

            h = this.x(2) - this.x(1);

            coll_points = zeros(2 * (n_y - 1), 1);
            for k = 1:(n_y - 1)
                map_to_coll = (1:2)' + 2 * (k - 1);
                coll_points(map_to_coll) = this.x(k) + h * ((1 / sqrt(3)) * [-1; 1] + 1) / 2;
            end

            nodes_to_coll_points = zeros(2 * (n_y - 1), n_y);
            for k = 1:(n_y - 1)
                map_to_coll = (1:2)' + 2 * (k - 1);
                nodes_to_coll_points(map_to_coll(1), k) = (coll_points(map_to_coll(1)) - this.x(k + 1)) / (this.x(k) - this.x(k + 1));
                nodes_to_coll_points(map_to_coll(1), k + 1) = (coll_points(map_to_coll(1)) - this.x(k)) / (this.x(k + 1) - this.x(k));
                nodes_to_coll_points(map_to_coll(2), k) = (coll_points(map_to_coll(2)) - this.x(k + 1)) / (this.x(k) - this.x(k + 1));
                nodes_to_coll_points(map_to_coll(2), k + 1) = (coll_points(map_to_coll(2)) - this.x(k)) / (this.x(k + 1) - this.x(k));
            end
            this.nodes_to_coll_points = nodes_to_coll_points;

            M = diag(4 * ones(1, n_y)) + diag(ones(1, n_y - 1), 1) + diag(ones(1, n_y - 1), -1);
            M(1, 1) = .5 * M(1, 1);
            M(end, end) = .5 * M(end, end);
            M = (1 / 6) * h * M;
            this.M = M;

            S = diag(2 * ones(1, n_y)) + (-1) * diag(ones(1, n_y - 1), 1) + (-1) * diag(ones(1, n_y - 1), -1);
            S(1, 1) = .5 * S(1, 1);
            S(end, end) = .5 * S(end, end);
            S = (1 / h) * S;
            this.S = S;

            this.dirichlet_bc = zeros(n_y, 1);
            this.dirichlet_bc(1) = 0;
            this.dirichlet_bc(n_y) = 0;

            this.forcing = this.x.^2;

            B = (10^3) * this.M;
            B(1, :) = 0 * B(1, :);
            B(end, :) = 0 * B(end, :);
            this.B = B;

        end

    end
end
