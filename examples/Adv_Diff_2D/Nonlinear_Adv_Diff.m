classdef Nonlinear_Adv_Diff < handle

    
    properties
        model;
        x;
        y;
    end
    
    methods
        function this = Nonlinear_Adv_Diff(adv_diff_obj)
            this.x = adv_diff_obj.pde_meshing.x;
            this.y = adv_diff_obj.pde_meshing.y;
            this.model = createpde;
            this.model.Geometry = adv_diff_obj.pde_meshing.model.Geometry;
            aCoef = @(region,state) adv_diff_obj.adv_coeff*(state.ux + state.uy); 
            specifyCoefficients(this.model,"m",0,"d",0,"c",adv_diff_obj.diff_coeff,"a",aCoef,"f",0);
            this.model.Mesh = adv_diff_obj.pde_meshing.model.Mesh;
            applyBoundaryCondition(this.model,"dirichlet","Edge",[3,4],"u",0);
            %this.model.BoundaryConditions = adv_diff_obj.pde_meshing.model.BoundaryConditions;
        end
        
        function [u] = State_Solve(this,z)
            control_fun = this.Map_z_to_Control_Fun(z);
            this.model.EquationCoefficients.CoefficientAssignments.f = @(region,~) control_fun(region.x,region.y);
            this.model.SolverOptions.MaxIterations = 1000;
            this.model.SolverOptions.ReportStatistics = 'on';
            this.model.SolverOptions.ResidualTolerance = 1.e-8;
            this.model.SolverOptions.MinStep = 1.e-9;
            result = solvepde(this.model);
            u = result.NodalSolution;
        end
 
        function [z] = Map_Control_Fun_to_z(this,control_fun)
            z=control_fun(this.x,this.y);
        end

        function [control_fun] = Map_z_to_Control_Fun(this,z)
           control_fun = scatteredInterpolant(this.x,this.y,z); 
        end
        
    end
end

