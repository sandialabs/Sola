classdef Diff_React_Hifi < handle

    
    properties
        model;
        x;
        y;
    end
    
    methods
        function obj = Diff_React_Hifi(pde_meshing)
            obj.x = pde_meshing.x;
            obj.y = pde_meshing.y;
            obj.model = createpde;
            obj.model.Geometry = pde_meshing.model.Geometry;
            aCoef = @(region,state) 5 + 0.05*(1+.1*sin(pi*region.x)).*state.u.^2; 
            c = 1;
            specifyCoefficients(obj.model,"m",0,"d",0,"c",c,"a",aCoef,"f",0);
            obj.model.BoundaryConditions = pde_meshing.model.BoundaryConditions;
            obj.model.Mesh = pde_meshing.model.Mesh;
        end
        
        function [u] = State_Solve(obj,control_fun)
            obj.model.EquationCoefficients.CoefficientAssignments.f = @(region,~) control_fun(region.x,region.y);
            obj.model.SolverOptions.MaxIterations = 1000;
            obj.model.SolverOptions.ReportStatistics = 'off';
            obj.model.SolverOptions.ResidualTolerance = 1.e-5;
            obj.model.SolverOptions.MinStep = 1.e-9;
            result = solvepde(obj.model);
            u = result.NodalSolution;
        end
 
        function [z] = Map_Control_Fun_to_z(obj,control_fun)
            z=control_fun(obj.x,obj.y);
        end

        function [control_fun] = Map_z_to_Control_Fun(obj,z)
           control_fun = scatteredInterpolant(obj.x,obj.y,z); 
        end
        
    end
end

