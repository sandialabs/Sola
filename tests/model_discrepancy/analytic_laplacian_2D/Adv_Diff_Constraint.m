classdef Adv_Diff_Constraint < Constraint

    properties
        adv_diff
        control_xlim
        control_ylim
        control_centers
        control_shape
        control_basis
        target_xlim
        target_ylim
        P_target
        m
        n
        M
        T
    end

    methods (Access = public)

        function [u] = State_Solve(this, z)

            rhs = this.adv_diff.M * this.control_basis * z;
            rhs(this.adv_diff.bnd_nodes) = 0;
            u = this.adv_diff.A \ rhs;
        end

        function [Mv] = c_u_Transpose_Inverse_Apply(this, v, u, z)
            Mv = (this.adv_diff.A') \ v;
        end

        function [Mv] = c_z_Transpose_Apply(this, v, u, z)
            A = -this.adv_diff.M * this.control_basis;
            A(this.adv_diff.bnd_nodes, :) = 0;
            Mv = A' * v;
        end

        function [Mv] = c_u_Inverse_Apply(this, v, u, z)
            Mv = this.adv_diff.A \ v;
        end

        function [Mv] = c_z_Apply(this, v, u, z)
            A = -this.adv_diff.M * this.control_basis;
            A(this.adv_diff.bnd_nodes, :) = 0;
            Mv = A * v;
        end

        function [Mv] = c_uu_Apply(this, v, u, z, lambda)
            Mv = zeros(this.m, 1);
        end

        function [Mv] = c_uz_Apply(this, v, u, z, lambda)
            Mv = zeros(this.m, 1);
        end

        function [Mv] = c_zu_Apply(this, v, u, z, lambda)
            Mv = zeros(this.n, 1);
        end

        function [Mv] = c_zz_Apply(this, v, u, z, lambda)
            Mv = zeros(this.n, 1);
        end

    end

    methods (Access = public)

        function this = Adv_Diff_Constraint(adv_diff)
            this = this@Constraint();
            this.adv_diff = adv_diff;
            this.m = size(this.adv_diff.A, 1);
            this.M = this.adv_diff.pde_meshing.M;
            x = this.adv_diff.pde_meshing.x;
            y = this.adv_diff.pde_meshing.y;

            this.control_xlim = [-0.8, 0.0];
            this.control_ylim = [-0.8, 0.0];
            this.target_xlim = [0.6, 0.7];
            this.target_ylim = [0.8, 0.9];

            I = find(x > this.target_xlim(1));
            I = intersect(I, find(x < this.target_xlim(2)));
            I = intersect(I, find(y > this.target_ylim(1)));
            I = intersect(I, find(y < this.target_ylim(2)));
            v = zeros(this.m, 1);
            v(I) = 1;
            this.P_target = diag(v);

            control_bandwidth = .2;
            x_control_nodes = this.control_xlim(1):control_bandwidth:this.control_xlim(2);
            y_control_nodes = this.control_ylim(1):control_bandwidth:this.control_xlim(2);
            nx = length(x_control_nodes);
            ny = length(y_control_nodes);

            this.n = nx * ny;
            this.control_centers = zeros(this.n, 2);
            this.control_centers(:, 1) = kron(ones(ny, 1), x_control_nodes');
            this.control_centers(:, 2) = kron(y_control_nodes', ones(nx, 1));
            this.control_shape = 30;

            this.control_basis = zeros(this.m, this.n);
            for k = 1:this.n
                tmp = [x, y] - this.control_centers(k, :);
                this.control_basis(:, k) = exp(-this.control_shape * sum(tmp.^2, 2));
            end

            this.T = 4 - 0 * x;

        end

        function [z_mesh] = Map_z_vec_to_mesh(this, z)
            z_mesh = this.control_basis * z;
        end

    end
end
