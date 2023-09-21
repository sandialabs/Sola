classdef Adv_Diff < Constrained_Optimization
    
    
    properties
        m;
        reg_coeff;
        reg_mat;
        T;
        diff_coeff;
        vel_coeff;
        robin_coeff;
        x;
        M;
        S;
        V;
        robin_bc;
    end
    
    methods (Access = public)
        
        %% Pure virtual functions for gradient computation
        
        function [val, grad_u, grad_z] = Objective(obj,u,z)
            val = (1/2)*(u-obj.T)'*obj.M*(u-obj.T) + (1/2)*(obj.reg_coeff)*z'*obj.reg_mat*z;
            grad_u = obj.M*(u-obj.T);
            grad_z = (obj.reg_coeff)*obj.reg_mat*z;
        end
        
        function [u] = State_Solve(obj,z)
            A = obj.diff_coeff*obj.S + obj.vel_coeff*obj.V + obj.robin_coeff*obj.robin_bc;
            b = (10^2)*obj.M*z;
            u = linsolve(A,b);
        end
        
        function [Mv] = c_u_Transpose_Inverse_Apply(obj,v,u,z)
            A = obj.diff_coeff*obj.S + obj.vel_coeff*obj.V + obj.robin_coeff*obj.robin_bc;
            Mv = linsolve(A',v);
        end
        
        function [Mv] = c_z_Transpose_Apply(obj,v,u,z)
            Mv = -(10^2)*obj.M'*v;
        end
        
        function [Mv] = c_u_Inverse_Apply(obj,v,u,z)
            A = obj.diff_coeff*obj.S + obj.vel_coeff*obj.V + obj.robin_coeff*obj.robin_bc;
            Mv = linsolve(A,v);
        end
        
        function [Mv] = c_z_Apply(obj,v,u,z)
            Mv = -(10^2)*obj.M*v;
        end
        
        function [Mv] = c_uu_Apply(obj,v,u,z,lambda)
            Mv = zeros(obj.m,1);
        end
        
        function [Mv] = c_uz_Apply(obj,v,u,z,lambda)
            Mv = zeros(obj.m,1);
        end
        
        function [Mv] = c_zu_Apply(obj,v,u,z,lambda)
            Mv = zeros(obj.m,1);
        end
        
        function [Mv] = c_zz_Apply(obj,v,u,z,lambda)
            Mv = zeros(obj.m,1);
        end
        
        function [Mv] = Objective_uu_Apply(obj,v,u,z)
            Mv = obj.M*v;
        end
        
        function [Mv] = Objective_uz_Apply(obj,v,u,z)
            Mv = zeros(obj.m,1);
        end
        
        function [Mv] = Objective_zu_Apply(obj,v,u,z)
            Mv = zeros(obj.m,1);
        end
        
        function [Mv] = Objective_zz_Apply(obj,v,u,z)
            Mv = obj.reg_coeff*obj.reg_mat*v;
        end
        
    end
    
    methods (Access = public)
        
        function obj = Adv_Diff(m,diff_coeff,vel_coeff,robin_coeff,reg_coeff)
            obj = obj@Constrained_Optimization();
            obj.m = m;
            obj.diff_coeff = diff_coeff;
            obj.vel_coeff = vel_coeff;
            obj.robin_coeff = robin_coeff;
            obj.x = linspace(0,1,m)';
            obj.reg_coeff = reg_coeff;
            obj.T = 50 - 30*(obj.x-0.5).^2;
            
            h = obj.x(2)-obj.x(1);
            
            M = diag(4*ones(1,m)) + diag(ones(1,m-1),1) + diag(ones(1,m-1),-1);
            M(1,1) = .5*M(1,1);
            M(end,end) = .5*M(end,end);
            M = (1/6)*h*M;
            obj.M = M;
            
            S = diag(2*ones(1,m)) + (-1)*diag(ones(1,m-1),1) + (-1)*diag(ones(1,m-1),-1);
            S(1,1) = .5*S(1,1);
            S(end,end) = .5*S(end,end);
            S = (1/h)*S;
            obj.S = S;
            
            V = diag(0*ones(1,m)) + (1/2)*diag(ones(1,m-1),1) + (-1/2)*diag(ones(1,m-1),-1);
            V(1,1) = -1/2;
            V(end,end) = 1/2;
            obj.V = V;
            
            robin_bc = zeros(m,m);
            robin_bc(1,1) = 1;
            robin_bc(end,end) = 1;
            obj.robin_bc = robin_bc;

            obj.reg_mat = M;
        end
        
    end
end

