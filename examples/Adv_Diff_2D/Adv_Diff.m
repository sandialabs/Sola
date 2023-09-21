classdef Adv_Diff < handle

    
    properties
        A;
        M;
        pde_meshing;
        diff_coeff;
        adv_coeff;
        bnd_nodes;
    end
    
    methods
        function obj = Adv_Diff(pde_meshing, diff_coeff, adv_coeff)
            model = createpde;
            model.Geometry = pde_meshing.model.Geometry;
            obj.diff_coeff = diff_coeff;
            obj.adv_coeff = adv_coeff;
            
            vel_field = @(x,y) adv_coeff*[ 1 + 0*x ; 1 + 0*y];
            fCoef = @(region,state) diag(vel_field(region.x,region.y)'*[state.ux; state.uy])';
            specifyCoefficients(model,"m",0,"d",0,"c",diff_coeff,"a",0,"f",fCoef);
            model.Mesh = pde_meshing.model.Mesh;

            m = size(model.Mesh.Nodes,2);
            state = struct;
            Adv_Op = zeros(m,m);
            for k = 1:m
                state.u = zeros(1,m);
                state.u(k) = 1;
                FEM = assembleFEMatrices(model,'domain',state);
                Adv_Op(:,k) = FEM.F;
            end
            obj.A = FEM.K + FEM.A + sparse(Adv_Op);
            obj.bnd_nodes = unique(union(find(pde_meshing.x==-1),find(pde_meshing.y==-1)));
            obj.A(obj.bnd_nodes,:) = 0;
            obj.A = obj.A + sparse(obj.bnd_nodes,obj.bnd_nodes,ones(length(obj.bnd_nodes),1),m,m);
            obj.M = pde_meshing.M;
            obj.pde_meshing = pde_meshing;
        end
        
        function [u] = State_Solve(obj,z)
           rhs = obj.M*z;
           rhs(obj.bnd_nodes) = 0; 
           u = obj.A\rhs; 
        end
        
        function [J_u] = State_Jacobian(obj)
           J_u = obj.A; 
        end
        
        function [J_z] = Control_Jacobian(obj)
            J_z = -obj.M;
            J_z(obj.bnd_nodes,:) = 0;
        end

    end
end

