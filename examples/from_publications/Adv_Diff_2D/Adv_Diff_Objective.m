%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Adv_Diff_Objective < Objective

    properties
        adv_diff
        reg_coeff
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

        function [val, grad_u, grad_z] = J(this, u, z)
            d = this.P_target * (u - this.T);
            val = (1 / 2) * d' * this.M * d + (1 / 2) * (this.reg_coeff) * (this.control_basis * z)' * this.M * (this.control_basis * z);
            grad_u = this.P_target' * this.M * d;
            grad_z = (this.reg_coeff) * this.control_basis' * this.M * (this.control_basis * z);
        end

        function [Mv] = J_uu_Apply(this, v, u, z)
            Mv = this.P_target' * this.M * this.P_target * v;
        end

        function [Mv] = J_uz_Apply(this, v, u, z)
            Mv = zeros(this.m, 1);
        end

        function [Mv] = J_zu_Apply(this, v, u, z)
            Mv = zeros(this.n, 1);
        end

        function [Mv] = J_zz_Apply(this, v, u, z)
            Mv = this.reg_coeff * this.control_basis' * this.M * this.control_basis * v;
        end

    end

    methods (Access = public)

        function this = Adv_Diff_Objective(adv_diff, reg_coeff)
            this = this@Objective();
            this.adv_diff = adv_diff;
            this.reg_coeff = reg_coeff;
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
