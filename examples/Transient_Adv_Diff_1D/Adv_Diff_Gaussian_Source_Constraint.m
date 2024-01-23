classdef Adv_Diff_Gaussian_Source_Constraint < Dynamic_Constraint

    properties
        A
        M
        B
        x
        num_space_control_nodes
        z_time_mesh
    end

    methods (Access = public)

        function [f, f_y, f_z] = f(this, y, z, t)
            % Extract the control for the given time.
            w = this.Temporal_Weights(t);
            zt = reshape(z, this.num_space_control_nodes, length(this.z_time_mesh)) * w;

            % Evaluate the RHS and its derivatives.
            f = linsolve(this.M, -this.A * y + this.M * this.B * zt);
            f_y = linsolve(this.M, -this.A);
            f_z = kron(w', this.B);
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

        function this = Adv_Diff_Gaussian_Source_Constraint(n_y, n, T, n_t, num_space_control_nodes)
            this = this@Dynamic_Constraint(n_y, n, T, n_t);

            % Peclet number
            Pe = 1;

            % Spatial domain
            this.x = linspace(0, 1, n_y)';
            dx = this.x(2) - this.x(1);

            % Mass matrix
            M = diag(4 * ones(1, n_y)) + diag(ones(1, n_y - 1), 1) + diag(ones(1, n_y - 1), -1);
            M(1, 1) = .5 * M(1, 1);
            M(end, end) = .5 * M(end, end);
            M = (1 / 6) * dx * M;

            % Stiffness matrix (diffusion)
            S = diag(2 * ones(1, n_y)) + (-1) * diag(ones(1, n_y - 1), 1) + (-1) * diag(ones(1, n_y - 1), -1);
            S(1, 1) = .5 * S(1, 1);
            S(end, end) = .5 * S(end, end);
            S = (1 / dx) * S;

            % Viscosity matrix (advection)
            V = diag(0 * ones(1, n_y)) + (1 / 2) * diag(ones(1, n_y - 1), 1) + (-1 / 2) * diag(ones(1, n_y - 1), -1);
            V(1, 1) = -1 / 2;
            V(end, end) = 1 / 2;

            % Spatially discretized PDE: M f(y,z) = -Ay + MBz,  A = S + Pe V
            A = S + Pe * V;
            this.A = A;
            this.M = M;

            % Control matrix (Gaussian)
            source_loc = linspace(0, 1, num_space_control_nodes)';
            B = zeros(n_y, num_space_control_nodes);
            for k = 1:num_space_control_nodes
                B(:, k) = exp(-200 * (this.x - source_loc(k)).^2);
            end
            B(1, :) = 0 * B(1, :);
            B(end, :) = 0 * B(end, :);

            % Store control properties
            this.B = B;
            this.num_space_control_nodes = num_space_control_nodes;
            this.z_time_mesh = linspace(0, T, n_t)';
            this.z_time_mesh = this.z_time_mesh(2:end);
        end

        function [w] = Temporal_Weights(this, t)
            w = (this.z_time_mesh - t) / (this.z_time_mesh(2) - this.z_time_mesh(1));
            Im = intersect(find(w <= 0), find(abs(w) <= 1));
            Ip = intersect(find(w > 0), find(abs(w) <= 1));
            I = find(abs(w) > 1);
            w(I) = 0;
            w(Im) = 1 + w(Im);
            w(Ip) = 1 - w(Ip);
        end

    end

end
