classdef Adv_Diff < handle

    properties
        A
        M
        pde_meshing
        diff_coeff
        adv_coeff
        bnd_nodes
    end

    methods

        function this = Adv_Diff(pde_meshing, diff_coeff, adv_coeff)
            model = createpde;
            model.Geometry = pde_meshing.model.Geometry;
            this.diff_coeff = diff_coeff;
            this.adv_coeff = adv_coeff;

            vel_field = @(x, y) adv_coeff * [1 + 0 * x; 1 + 0 * y];
            fCoef = @(region, state) diag(vel_field(region.x, region.y)' * [state.ux; state.uy])';
            specifyCoefficients(model, "m", 0, "d", 0, "c", diff_coeff, "a", 0, "f", fCoef);
            model.Mesh = pde_meshing.model.Mesh;

            m = size(model.Mesh.Nodes, 2);
            state = struct;
            Adv_Op = zeros(m, m);
            for k = 1:m
                state.u = zeros(1, m);
                state.u(k) = 1;
                FEM = assembleFEMatrices(model, 'domain', state);
                Adv_Op(:, k) = FEM.F;
            end
            this.A = FEM.K + FEM.A + sparse(Adv_Op);
            this.bnd_nodes = unique(find(pde_meshing.y == -1));
            this.A(this.bnd_nodes, :) = 0;
            this.A = this.A + sparse(this.bnd_nodes, this.bnd_nodes, ones(length(this.bnd_nodes), 1), m, m);
            this.M = pde_meshing.M;
            this.pde_meshing = pde_meshing;
        end

        function [u] = State_Solve(this, z)
            rhs = this.M * z;
            rhs(this.bnd_nodes) = 0;
            u = this.A \ rhs;
        end

        function [J_u] = State_Jacobian(this)
            J_u = this.A;
        end

        function [J_z] = Control_Jacobian(this)
            J_z = -this.M;
            J_z(this.bnd_nodes, :) = 0;
        end

    end
end
