classdef Linear_1D_Finite_Elements < handle

    properties
        m
        x
        M
        S
        nodes_to_coll_points
        coll_point_integration
    end

    methods

        function this = Linear_1D_Finite_Elements(m)
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

            coll_point_integration = zeros(m, 2 * (m - 1));
            a = (h / 2) * ((1 - 1 / sqrt(3)) / 2);
            A = (h / 2) * ((1 + 1 / sqrt(3)) / 2);
            coll_point_integration(1, 1:2) = [A; a];
            coll_point_integration(end, end - 1:end) = [a; A];
            for k = 2:(m - 1)
                coll_point_integration(k, (2 * (k - 2) + 1):(2 * k)) = [a; A; A; a];
            end
            this.coll_point_integration = coll_point_integration;

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

        end

    end
end
