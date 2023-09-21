classdef PDE_Meshing < handle

    properties
        S; % Stiffness matrix
        M; % Mass matrix
        x; % x nodes
        y; % y nodes
        model; % Model object which contains mesh and geometry
    end
    
    methods
        function obj = PDE_Meshing(Hmax)
            
            model = createpde;
            
            %geometryFromEdges(model,@circleg);
            geometryFromEdges(model,@squareg);
            if false
                pdegplot(model,"EdgeLabels","on")
                axis equal
                title("Geometry with Edge Labels")
            end
            
            generateMesh(model,"Hmax",Hmax);
            if false
                figure;
                pdemesh(model);
                axis equal
            end
            
            specifyCoefficients(model,"m",0,"d",0,"c",1,"a",1,"f",0);
            FEM = assembleFEMatrices(model,'KA');
            
            obj.S = FEM.K;
            obj.M = FEM.A;
            obj.x = model.Mesh.Nodes(1,:)';
            obj.y = model.Mesh.Nodes(2,:)';
            obj.model = model;
            
        end
        
        function [] = Plot_Field(obj,u,name)
            xl = linspace(min(obj.x),max(obj.x),100)';
            yl = linspace(min(obj.y),max(obj.y),100)';
            [X,Y] = meshgrid(xl,yl);
            F = scatteredInterpolant(obj.x,obj.y,u);
            f = F(X,Y);
            figure,
            surf(X,Y,f)
            view(2)
            shading interp
            colorbar()
            title(name)
        end
        
    end
end

