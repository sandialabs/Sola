classdef Example_1 < Constrained_Optimization 
    
    %% Problem description
    % $$\min_{z \in \mathbf{R}^2} J(u,z) = || S(z) - (7,1,4)^T ||^2 + (z_1-8)^2 + (z_2-8)^2 + (u(1)z(1)-56)^2$$
    %
    % s.t 
    %
    % $$u = S(z)$$ 
    %
    % solves
    %
    % $$u_1 + u_2     + 0     = z_1$$
    %
    % $$0   + z_1u_2 + 0     = z_2$$
    %
    % $$0   +  0      + u_3^3 = z_2^2$$
    %
    % We derive off <Constrained_Optimization.html Constrained_Optimization>
    
    methods (Access = public)

        function obj = Example_1( )

        end
        
        %% Instantiation of base class pure virtual functions for gradient computation
        function [val,grad_u, grad_z] = Objective(obj,u,z)
            val = norm(u-[7;1;4])^2 + (z(1)-8).^2 + (z(2)-8).^2 + (u(1)*z(1)-56).^2;
            grad_u = 2*(u-[7;1;4]);
            grad_u(1) = grad_u(1) + 2*(u(1)*z(1)-56)*z(1);
            grad_z = zeros(2,1);
            grad_z(1) = 2*(z(1)-8) + 2*(u(1)*z(1)-56)*u(1);
            grad_z(2) = 2*(z(2)-8);
        end
        
        function [u] = State_Solve(obj,z)
            u = zeros(3,1);
            u(1) = z(1) - z(2)/z(1);
            u(2) = z(2)/z(1);
            u(3) = z(2).^(2/3);
        end
        
        function [Mv] = c_u_Transpose_Inverse_Apply(obj,v,u,z) 
            A = [1 1 0; 0 z(1) 0; 0 0 3*u(3).^2];
            Mv = linsolve(A',v);
        end

        function [Mv] = c_z_Transpose_Apply(obj,v,u,z)
            A = [ -1 , u(2) , 0 ; 0 , -1 , -2*z(2) ];
            Mv = A*v;
        end
        
        %% Instantiation of base class pure virtual functions for hessian-vector product computation
        function [Mv] = c_u_Inverse_Apply(obj,v,u,z) 
           A = [1 , 1 , 0 ; 0 , z(1) , 0 ; 0 , 0 , 3*u(3)^2 ];
           Mv = linsolve(A,v);
        end
        
        function [Mv] = c_z_Apply(obj,v,u,z) 
            A = [ -1 , u(2) , 0 ; 0 , -1 , -2*z(2) ];
            Mv = A'*v;
        end
        
        function [Mv] = c_uu_Apply(obj,v,u,z,lambda) 
            A = zeros(3,3);
            A(3,3) = 6*lambda(3)*u(3);
            Mv = A*v;
        end
        
        function [Mv] = c_uz_Apply(obj,v,u,z,lambda) 
            A = zeros(length(u),length(z));
            A(2,1) = lambda(2);
            Mv = A*v;
        end
        
        function [Mv] = c_zu_Apply(obj,v,u,z,lambda) 
            A = zeros(length(z),length(u));
            A(1,2) = lambda(2);
            Mv = A*v;
        end
        
        function [Mv] = c_zz_Apply(obj,v,u,z,lambda)
            A = zeros(2,2);
            A(2,2) = -2*lambda(3);
            Mv = A*v;
        end
        
        function [Mv] = Objective_uu_Apply(obj,v,u,z)
            A = 2*eye(3);
            A(1,1) = A(1,1) + 2*z(1)^2;
            Mv = A*v;
        end
        
        function [Mv] = Objective_uz_Apply(obj,v,u,z) 
            A = zeros(length(u),length(z));
            A(1,1) = 2*(u(1)*z(1)-56) + 2*u(1)*z(1);
            Mv = A*v;
        end
        
        function [Mv] = Objective_zu_Apply(obj,v,u,z) 
            A = zeros(length(z),length(u));
            A(1,1) = 2*(u(1)*z(1)-56) + 2*u(1)*z(1);
            Mv = A*v;
        end
        
        function [Mv] = Objective_zz_Apply(obj,v,u,z) 
            A = 2*eye(2);
            A(1,1) = A(1,1) + 2*u(1)^2;
            Mv = A*v;
        end
        
    end

end

