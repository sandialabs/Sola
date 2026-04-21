%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Diff_React_Lofi < handle

    properties
        A
        M
        S
        pde_meshing
    end

    methods

        function this = Diff_React_Lofi(pde_meshing)
            model = createpde;
            model.Geometry = pde_meshing.model.Geometry;

            aCoef = @(region, ~) 5 * 1 + 0 * region.x;
            c = 1;
            specifyCoefficients(model, "m", 0, "d", 0, "c", c, "a", aCoef, "f", 0);
            model.BoundaryConditions = pde_meshing.model.BoundaryConditions;
            model.Mesh = pde_meshing.model.Mesh;
            FEM = assembleFEMatrices(model, 'domain');
            this.A = FEM.K + FEM.A;
            this.M = pde_meshing.M;
            this.S = FEM.K;
            this.pde_meshing = pde_meshing;
        end

        function [u] = State_Solve(this, z)
            u = this.A \ (this.M * z);
        end

        function [J_u] = State_Jacobian(this)
            J_u = this.A;
        end

        function [J_z] = Control_Jacobian(this)
            J_z = -this.M;
        end

    end
end
