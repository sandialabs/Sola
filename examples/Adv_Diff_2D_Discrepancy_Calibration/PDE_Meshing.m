classdef PDE_Meshing < handle

    properties
        S  % Stiffness matrix
        M  % Mass matrix
        x  % x nodes
        y  % y nodes
        model  % Model object which contains mesh and geometry
    end

    methods

        function this = PDE_Meshing(Hmax)

            model = createpde;

            R1 = [3,4,0,1,1,0,1,1,0,0]';
            g = decsg(R1);
            geometryFromEdges(model,g);
            %geometryFromEdges(model, @squareg);
            if false
                pdegplot(model, "EdgeLabels", "on");
                axis equal;
                title("Geometry with Edge Labels");
            end

            generateMesh(model, "Hmax", Hmax);
            if false
                figure;
                pdemesh(model);
                axis equal;
            end

            specifyCoefficients(model, "m", 0, "d", 0, "c", 1, "a", 1, "f", 0);
            FEM = assembleFEMatrices(model, 'KA');

            this.S = FEM.K;
            this.M = FEM.A;
            this.x = model.Mesh.Nodes(1, :)';
            this.y = model.Mesh.Nodes(2, :)';
            this.model = model;

        end

        function [] = Plot_Field(this, u, name)
            xl = linspace(min(this.x), max(this.x), 100)';
            yl = linspace(min(this.y), max(this.y), 100)';
            [X, Y] = meshgrid(xl, yl);
            F = scatteredInterpolant(this.x, this.y, u);
            f = F(X, Y);
            figure;
            surf(X, Y, f);
            view(2);
            shading interp;
            colorbar();
            title(name);
        end

    end
end
