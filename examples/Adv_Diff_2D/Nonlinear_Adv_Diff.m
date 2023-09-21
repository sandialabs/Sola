classdef Nonlinear_Adv_Diff < handle

    
    properties
        model;
        x;
        y;
    end
    
    methods
        function obj = Nonlinear_Adv_Diff(adv_diff_obj)
            obj.x = adv_diff_obj.pde_meshing.x;
            obj.y = adv_diff_obj.pde_meshing.y;
            obj.model = createpde;
            obj.model.Geometry = adv_diff_obj.pde_meshing.model.Geometry;
            aCoef = @(region,state) adv_diff_obj.adv_coeff*(state.ux + state.uy); 
            specifyCoefficients(obj.model,"m",0,"d",0,"c",adv_diff_obj.diff_coeff,"a",aCoef,"f",0);
            obj.model.Mesh = adv_diff_obj.pde_meshing.model.Mesh;
            applyBoundaryCondition(obj.model,"dirichlet","Edge",[3,4],"u",0);
            %obj.model.BoundaryConditions = adv_diff_obj.pde_meshing.model.BoundaryConditions;
        end
        
        function [u] = State_Solve(obj,z)
            control_fun = obj.Map_z_to_Control_Fun(z);
            obj.model.EquationCoefficients.CoefficientAssignments.f = @(region,~) control_fun(region.x,region.y);
            obj.model.SolverOptions.MaxIterations = 1000;
            obj.model.SolverOptions.ReportStatistics = 'on';
            obj.model.SolverOptions.ResidualTolerance = 1.e-8;
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

