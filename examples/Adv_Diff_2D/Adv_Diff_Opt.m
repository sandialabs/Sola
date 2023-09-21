classdef Adv_Diff_Opt < Constrained_Optimization
    
    
    properties
        adv_diff;
        reg_coeff;
        control_xlim;
        control_ylim;
        control_centers;
        control_shape;
        control_basis;
        target_xlim;
        target_ylim;
        P_target;
        m;
        n;
        M;
        T;
    end
    
    methods (Access = public)
        
        function [val, grad_u, grad_z] = Objective(obj,u,z)
            d = obj.P_target*(u-obj.T);
            val = (1/2)*d'*obj.M*d + (1/2)*(obj.reg_coeff)*(obj.control_basis*z)'*obj.M*(obj.control_basis*z);
            grad_u = obj.P_target'*obj.M*d;
            grad_z = (obj.reg_coeff)*obj.control_basis'*obj.M*(obj.control_basis*z);
        end
        
        function [u] = State_Solve(obj,z)
            
            rhs = obj.adv_diff.M*obj.control_basis*z;
            rhs(obj.adv_diff.bnd_nodes) = 0;
            u = obj.adv_diff.A\rhs;
        end
        
        function [Mv] = c_u_Transpose_Inverse_Apply(obj,v,u,z)
            Mv = (obj.adv_diff.A')\v;
        end
        
        function [Mv] = c_z_Transpose_Apply(obj,v,u,z)
            A = -obj.adv_diff.M*obj.control_basis;
            A(obj.adv_diff.bnd_nodes,:) = 0;
            Mv = A'*v;
        end
        
        function [Mv] = c_u_Inverse_Apply(obj,v,u,z)
            Mv = obj.adv_diff.A\v;
        end
        
        function [Mv] = c_z_Apply(obj,v,u,z)
            A = -obj.adv_diff.M*obj.control_basis;
            A(obj.adv_diff.bnd_nodes,:) = 0;
            Mv = A*v;
        end
        
        function [Mv] = c_uu_Apply(obj,v,u,z,lambda)
            Mv = zeros(obj.m,1);
        end
        
        function [Mv] = c_uz_Apply(obj,v,u,z,lambda)
            Mv = zeros(obj.m,1);
        end
        
        function [Mv] = c_zu_Apply(obj,v,u,z,lambda)
            Mv = zeros(obj.n,1);
        end
        
        function [Mv] = c_zz_Apply(obj,v,u,z,lambda)
            Mv = zeros(obj.n,1);
        end
        
        function [Mv] = Objective_uu_Apply(obj,v,u,z)
            Mv = obj.P_target'*obj.M*obj.P_target*v;
        end
        
        function [Mv] = Objective_uz_Apply(obj,v,u,z)
            Mv = zeros(obj.m,1);
        end
        
        function [Mv] = Objective_zu_Apply(obj,v,u,z)
            Mv = zeros(obj.n,1);
        end
        
        function [Mv] = Objective_zz_Apply(obj,v,u,z)
            Mv = obj.reg_coeff*obj.control_basis'*obj.M*obj.control_basis*v;
        end
        
    end
    
    methods (Access = public)
         
        function obj = Adv_Diff_Opt(adv_diff,reg_coeff)
            obj = obj@Constrained_Optimization();
            obj.adv_diff = adv_diff;
            obj.reg_coeff = reg_coeff;
            obj.m = size(obj.adv_diff.A,1);
            obj.M = obj.adv_diff.pde_meshing.M;
            x = obj.adv_diff.pde_meshing.x;
            y = obj.adv_diff.pde_meshing.y;
            
            obj.control_xlim = [-0.8,0.0];
            obj.control_ylim = [-0.8,0.0];
            obj.target_xlim = [0.6,0.7];
            obj.target_ylim = [0.8,0.9];
            
            I = find(x > obj.target_xlim(1));
            I = intersect(I,find(x < obj.target_xlim(2)));
            I = intersect(I,find(y > obj.target_ylim(1)));
            I = intersect(I,find(y < obj.target_ylim(2)));
            v = zeros(obj.m,1);
            v(I) = 1;
            obj.P_target = diag(v);
            
            control_bandwidth = .2;
            x_control_nodes = obj.control_xlim(1):control_bandwidth:obj.control_xlim(2);
            y_control_nodes = obj.control_ylim(1):control_bandwidth:obj.control_xlim(2);
            nx = length(x_control_nodes);
            ny = length(y_control_nodes);
            
            obj.n = nx*ny;
            obj.control_centers = zeros(obj.n,2);
            obj.control_centers(:,1) = kron(ones(ny,1),x_control_nodes');
            obj.control_centers(:,2) = kron(y_control_nodes',ones(nx,1));
            obj.control_shape = 30;
            
            obj.control_basis = zeros(obj.m,obj.n);
            for k = 1:obj.n
                tmp = [x,y] - obj.control_centers(k,:);
                obj.control_basis(:,k) = exp(-obj.control_shape*sum(tmp.^2,2));
            end
            
            obj.T = 4 - 0*x;
            
        end
        
        function [z_mesh] = Map_z_vec_to_mesh(obj,z)
            z_mesh = obj.control_basis*z;
        end
        
    end
end

