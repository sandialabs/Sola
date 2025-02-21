classdef Diff_Constraint < Constraint

    properties
        m
        M
        S
        A
        x
        y
        bnd_nodes
        pde_meshing
    end

    methods (Access = public)

        function [u] = State_Solve(this, z)

            rhs = this.M * z;
            rhs(this.bnd_nodes) = 0;
            u = this.A \ rhs;
        end

        function [Mv] = c_u_Transpose_Inverse_Apply(this, v, u, z)
            Mv = (this.A') \ v;
        end

        function [Mv] = c_z_Transpose_Apply(this, v, u, z)
            B = -this.M;
            B(this.bnd_nodes, :) = 0;
            Mv = B' * v;
        end

        function [Mv] = c_u_Inverse_Apply(this, v, u, z)
            Mv = this.A \ v;
        end

        function [Mv] = c_z_Apply(this, v, u, z)
            B = -this.M;
            B(this.bnd_nodes, :) = 0;
            Mv = B * v;
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

        function this = Diff_Constraint(pde_meshing, diff_coeff)
            this = this@Constraint();

            this.M = pde_meshing.M;
            this.S = pde_meshing.S;
            this.m = size(this.M, 1);
            this.x = pde_meshing.x;
            this.y = pde_meshing.y;
            this.bnd_nodes = unique(find(this.y == -1));

            this.A = diff_coeff * this.S;
            this.A(this.bnd_nodes, :) = 0;
            this.A = this.A + sparse(this.bnd_nodes, this.bnd_nodes, ones(length(this.bnd_nodes), 1), this.m, this.m);
            this.pde_meshing = pde_meshing;
        end

    end
end
