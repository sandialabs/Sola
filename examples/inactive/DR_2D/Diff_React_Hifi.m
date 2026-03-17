classdef Diff_React_Hifi < handle

    properties
        model
        x
        y
    end

    methods

        function this = Diff_React_Hifi(pde_meshing)
            this.x = pde_meshing.x;
            this.y = pde_meshing.y;
            this.model = createpde;
            this.model.Geometry = pde_meshing.model.Geometry;
            aCoef = @(region, state) 5 + 0.1 * (1 + .1 * sin(pi * region.x)) .* state.u.^2;
            c = 1;
            specifyCoefficients(this.model, "m", 0, "d", 0, "c", c, "a", aCoef, "f", 0);
            this.model.BoundaryConditions = pde_meshing.model.BoundaryConditions;
            this.model.Mesh = pde_meshing.model.Mesh;
        end

        function [u] = State_Solve(this, control_fun)
            this.model.EquationCoefficients.CoefficientAssignments.f = @(region, ~) control_fun(region.x, region.y);
            this.model.SolverOptions.MaxIterations = 1000;
            this.model.SolverOptions.ReportStatistics = 'off';
            this.model.SolverOptions.ResidualTolerance = 1.e-5;
            this.model.SolverOptions.MinStep = 1.e-9;
            result = solvepde(this.model);
            u = result.NodalSolution;
        end

        function [z] = Map_Control_Fun_to_z(this, control_fun)
            z = control_fun(this.x, this.y);
        end

        function [control_fun] = Map_z_to_Control_Fun(this, z)
            control_fun = scatteredInterpolant(this.x, this.y, z);
        end

    end
end
