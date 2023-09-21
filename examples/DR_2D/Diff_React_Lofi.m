classdef Diff_React_Lofi < handle

    
    properties
        A;
        M;
        pde_meshing;
    end
    
    methods
        function obj = Diff_React_Lofi(pde_meshing)
            model = createpde;
            model.Geometry = pde_meshing.model.Geometry;
            
            aCoef = @(region,~) 5*1 + 0*region.x;
            c = 1;
            specifyCoefficients(model,"m",0,"d",0,"c",c,"a",aCoef,"f",0);
            model.BoundaryConditions = pde_meshing.model.BoundaryConditions;
            model.Mesh = pde_meshing.model.Mesh;
            FEM = assembleFEMatrices(model,'domain');
            obj.A = FEM.K+FEM.A;
            obj.M = pde_meshing.M;
            obj.pde_meshing = pde_meshing;
        end
        
        function [u] = State_Solve(obj,z)
           u = obj.A\(obj.M*z); 
        end
        
        function [J_u] = State_Jacobian(obj)
           J_u = obj.A; 
        end
        
        function [J_z] = Control_Jacobian(obj)
            J_z = -obj.M;
        end

    end
end

