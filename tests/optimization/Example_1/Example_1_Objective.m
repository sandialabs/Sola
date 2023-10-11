classdef Example_1_Objective < Objective 
    
    %% Problem description
    % \min_{z \in R^2} J(u,z) = || S(z) - (7,1,4)^T ||^2 + (z_1-8)^2 + (z_2-8)^2 + (u(1)z(1)-56)^2
    %
    % s.t 
    %
    % u = S(z)
    %
    % solves
    %
    % u_1 + u_2     + 0     = z_1
    %
    % 0   + z_1u_2 + 0     = z_2
    %
    % 0   +  0      + u_3^3 = z_2^2
    
    methods (Access = public)

        function this = Example_1_Objective( )

        end
        
        function [val,grad_u, grad_z] = J(this,u,z)
            val = norm(u-[7;1;4])^2 + (z(1)-8).^2 + (z(2)-8).^2 + (u(1)*z(1)-56).^2;
            grad_u = 2*(u-[7;1;4]);
            grad_u(1) = grad_u(1) + 2*(u(1)*z(1)-56)*z(1);
            grad_z = zeros(2,1);
            grad_z(1) = 2*(z(1)-8) + 2*(u(1)*z(1)-56)*u(1);
            grad_z(2) = 2*(z(2)-8);
        end
        
        function [Mv] = J_uu_Apply(this,v,u,z)
            A = 2*eye(3);
            A(1,1) = A(1,1) + 2*z(1)^2;
            Mv = A*v;
        end
        
        function [Mv] = J_uz_Apply(this,v,u,z) 
            A = zeros(length(u),length(z));
            A(1,1) = 2*(u(1)*z(1)-56) + 2*u(1)*z(1);
            Mv = A*v;
        end
        
        function [Mv] = J_zu_Apply(this,v,u,z) 
            A = zeros(length(z),length(u));
            A(1,1) = 2*(u(1)*z(1)-56) + 2*u(1)*z(1);
            Mv = A*v;
        end
        
        function [Mv] = J_zz_Apply(this,v,u,z) 
            A = 2*eye(2);
            A(1,1) = A(1,1) + 2*u(1)^2;
            Mv = A*v;
        end
        
    end

end

