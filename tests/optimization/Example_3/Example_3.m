classdef Example_3 < Constrained_ODE_Optimization
    
    % Solve the optimiation problem
    % min_{z} J(S(z),z) = int_0^T g(S(z)(t))dt + R(z)
    % where S(z) solves the odinary differential equation
    % du/dt = [u3/u2, 2*u1, 3*u1*u2]
    % u(0) = [z_1^2 , z_2^3, z_3^4]
    % g(u) = (u_1-exp(t))^2 + (u_2-exp(2*t))^2 + (u_3-exp(3*t))^2
    % R(z) = (z_1-1)^2 + (z_2-1)^2 + (z_3-1)^2
    
    properties

    end
    
    methods (Access = public)
        
         %% Instantiation of base class pure virtual functions for gradient computation
        function [val, grad_u] = Time_Instance_Objective(obj,u,t) 
            val = (u(1)-exp(t))^2 + (u(2)-exp(2*t))^2 + (u(3)-exp(3*t))^2;
            grad_u = zeros(3,1);
            grad_u(1) = 2*(u(1)-exp(t));
            grad_u(2) = 2*(u(2)-exp(2*t));
            grad_u(3) = 2*(u(3)-exp(3*t));
        end
        
        function [val,grad_z] = Regularization_Objective(obj,z) 
            val = (z(1)-1)^2 + (z(2)-1)^2 + (z(3)-1)^2;
            grad_z = 2*(z-1);
        end
        
        function [f, f_u, f_z] = Time_Instance_RHS(obj,u,z,t) 
            f = [ u(3)/u(2) ; 2*u(1) ; 3*u(1)*u(2) ];
            f_u = [ 0 , -u(3)/u(2)^2 , 1/u(2) ; 2 , 0 , 0 ; 3*u(2) , 3*u(1) , 0 ];
            f_z = zeros(3,3);
        end
        
        function [h, h_z] = Initial_Condition(obj,z)
            h = [z(1)^2 ; z(2)^3 ; z(3)^4];
            h_z = [ 2*z(1) , 0 , 0 ; 0 , 3*z(2).^2 , 0 ; 0 , 0 , 4*z(3)^3 ];
        end
        
        %% Instantiation of base class pure virtual functions for hessian-vector product computation
        function [Mv] = Time_Instance_Objective_uu_Apply(obj,v,u,t) 
            A = 2*eye(3);
            Mv = A*v;
        end
        
        function [Mv] = Regularization_Objective_zz_Apply(obj,v,z) 
            A = 2*eye(3);
            Mv = A*v;
        end
        
        function [Mv] = Time_Instance_RHS_uu_Apply(obj,v,u,z,t,lambda) 
            A = [ 0 , 3*lambda(3) , 0 ; 3*lambda(3) , 2*lambda(1)*u(3)/u(2)^3 , -lambda(1)/u(2)^2 ; 0 , -lambda(1)/u(2)^2 , 0];
            Mv = A*v;
        end
        
        function [Mv] = Time_Instance_RHS_uz_Apply(obj,v,u,z,t,lambda) 
            num_vecs = size(v,2);
            Mv = zeros(obj.m,num_vecs);
        end
        
        function [Mv] = Time_Instance_RHS_zu_Apply(obj,v,u,z,t,lambda) 
            num_vecs = size(v,2);
            Mv = zeros(length(z),num_vecs);
        end
        
        function [Mv] = Time_Instance_RHS_zz_Apply(obj,v,u,z,t,lambda) 
            num_vecs = size(v,2);
            Mv = zeros(length(z),num_vecs);
        end
        
        function [Mv] = Initial_Condition_zz_Apply(obj,v,z,lambda) 
            A = diag([2*lambda(1) ; 6*lambda(2)*z(2) ; 12*lambda(3)*z(3)^2 ]);
            Mv = A*v;
        end
        
    end
    
    methods (Access = public)
        function obj = Example_3(m,n,T,N)
            obj = obj@Constrained_ODE_Optimization(m,n,T,N);
        end
       
    end
      
end

