 classdef Mass_Spring_Coupled < Constrained_ODE_Optimization
    
    properties
        m1;
        m2;
        k1;
        k2;
        k3;
        f2;
        reg_coeff;
        P_z;
    end
    
    methods 
        
        function obj = Mass_Spring_Coupled(T,N)
            m = 4;
            n = 1;
            obj = obj@Constrained_ODE_Optimization(m,n,T,N);
            
            obj.m1 = 1;
            obj.m2 = 10;
            obj.k1 = 1;
            obj.k2 = 1;
            obj.k3 = 1;
            obj.f2 = 0;
            obj.reg_coeff = 1.e-6;
            
            P_z = eye(N);
            P_z = P_z(:,2:end);
            P_z(1,1) = 1;
            obj.P_z = P_z;
        end
        
        function [w] = Temporal_Weights(obj,t)
            w = (obj.t_mesh-t)/(obj.t_mesh(2)-obj.t_mesh(1));
            Im = intersect(find(w<=0),find(abs(w)<=1));
            Ip = intersect(find(w>0),find(abs(w)<=1));
            w(abs(w)>1) = 0;
            w(Im) = 1+w(Im);
            w(Ip) = 1-w(Ip);
        end
        
        function [val] = target(obj,t)
            val = 5*t.^2;
        end
        
        function [val, grad_u] = Time_Instance_Objective(obj,u,t)
            val = 0.5*(u(1) - obj.target(t)).^2;
            grad_u = zeros(4,1);
            grad_u(1) = u(1)-obj.target(t);
        end
        
        function [val,grad_z] = Regularization_Objective(obj,z)
            val = 0.5*obj.reg_coeff*(obj.w'*(obj.P_z*z).^2);
            grad_z = obj.reg_coeff*obj.P_z'*(obj.w.*(obj.P_z*z));
        end
        
        % ODE system with four states u=(x_1,v_1,x_2,v_2)
        % x_1' = v_1
        % v_1' = (1/m_1)*( k_2*x_2 - (k_1+k_2)*x_1 + f_1(z) )
        % x_2' = v_2
        % v_2' = (1/m_2)*( k_2*x_1 - (k_2+k_3)*x_2 + f_2 )
        function [f, f_u, f_z] = Time_Instance_RHS(obj,u,z,t)
            x1 = u(1);
            v1 = u(2);
            x2 = u(3);
            v2 = u(4);

            coeffs = obj.Temporal_Weights(t);
            f1 = (obj.P_z*z)'*coeffs;
            
            f = zeros(4,1);
            f(1) = v1;
            f(2) = (1/obj.m1)*( obj.k2*x2 - (obj.k1+obj.k2)*x1 + f1 );
            f(3) = v2;
            f(4) = (1/obj.m2)*( obj.k2*x1 - (obj.k2+obj.k3)*x2 + obj.f2);
            
            f_u = zeros(4,4);
            f_u(1,2) = 1;
            f_u(2,1) = -(1/obj.m1)*(obj.k1+obj.k2);
            f_u(2,3) = (1/obj.m1)*obj.k2;
            f_u(3,4) = 1;
            f_u(4,1) = (1/obj.m2)*obj.k2;
            f_u(4,3) = -(1/obj.m2)*(obj.k2+obj.k3);
            
            f_z = zeros(4,size(obj.P_z,2));
            f_z(2,:) = (1/obj.m1)*obj.P_z'*coeffs;
        end
        
        function [h, h_z] = Initial_Condition(obj,z)
            h = zeros(4,1);
            h_z = zeros(4,size(obj.P_z,2));
        end
        
        function [Mv] = Time_Instance_Objective_uu_Apply(obj,v,u,t)
            Mv = zeros(size(v));
            Mv(1,:) = v(1,:);
        end
        
        function [Mv] = Regularization_Objective_zz_Apply(obj,v,z)
            Mv = obj.reg_coeff*obj.P_z'*diag(obj.w)*obj.P_z*v;
        end
        
        function [Mv] = Time_Instance_RHS_uu_Apply(obj,v,u,z,t,lambda)
            Mv = 0*v;
        end
        
        function [Mv] = Time_Instance_RHS_uz_Apply(obj,v,u,z,t,lambda)
            Mv = zeros(4,size(v,2));
        end
        
        function [Mv] = Time_Instance_RHS_zu_Apply(obj,v,u,z,t,lambda)
            Mv = zeros(length(z),size(v,2));
        end
        
        function [Mv] = Time_Instance_RHS_zz_Apply(obj,v,u,z,t,lambda)
            Mv = 0*v;
        end
        
        function [Mv] = Initial_Condition_zz_Apply(obj,v,z,lambda)
            Mv = 0*v;
        end
        
    end
    
    
    
end


