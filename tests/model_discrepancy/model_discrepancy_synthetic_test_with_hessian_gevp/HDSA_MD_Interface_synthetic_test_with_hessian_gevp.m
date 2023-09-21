classdef HDSA_MD_Interface_synthetic_test_with_hessian_gevp < HDSA_MD_Interface
    
    properties
        m; % Mesh resolution
        x; % Mesh nodes on [0,1]
        S; % Stiffness matrix
        M; % Mass matrix
        W_d; % Discrepancy precision matrix
        W_u; % State weighting matrix
        W_z; % Control weighting matrix
    end
    
    methods (Access = public)

        function [u_out] = Apply_W_d(obj,u_in)
            u_out = obj.W_d*u_in;
        end
        
        function [z_out] = Apply_W_z_Inverse(obj,z_in)
           z_out = linsolve(obj.W_z,z_in); 
        end
        
        function [z_out] = Apply_W_z_Inverse_Factor(obj,z_in)
            R = chol(obj.W_z);
            z_out = linsolve(R,z_in);
        end
        
        function [z_out] = Apply_W_z(obj,z_in)
           z_out = obj.W_z*z_in; 
        end
        
        function [u_out] = Apply_W_u_Plus_scalar_W_d_Inverse(obj,u_in,scalar)
            u_out = linsolve(obj.W_u + scalar*obj.W_d,u_in);
        end
        
        function [u_out] = Apply_W_u_Plus_scalar_W_d_Inverse_Factor(obj,u_in,scalar)
           R = chol(obj.W_u+scalar*obj.W_d);
           u_out = linsolve(R,u_in);
        end
        
        function [u_out] = Apply_W_u_Inverse(obj,u_in)
            u_out = linsolve(obj.W_u,u_in);
        end
        
        function [u_out] = Apply_W_u_Inverse_Factor(obj,u_in)
            R = chol(obj.W_u);
            u_out = linsolve(R,u_in);
        end
        
        function [z_out] = Apply_Solution_Operator_z_Jacobian_Transpose(obj,u_in,z)
            z_out = 3*diag(z.^2)*u_in;
        end
            
        % This implementation assumes that it is evaluated at the optimal z so that 
        % the adjoint=0, a more general impl a term multiplied by the adjoint variable
        function [z_out] = Apply_RS_Hessian(obj,z_in,z)
            tmp1 = 3*diag(z.^2)*z_in;
            tmp2 = obj.M*tmp1;
            z_out = 3*diag(z.^2)*tmp2;
        end
        
        function [grad_u] = Misfit_Gradient(obj,u,z)
            grad_u = obj.M*(u-(1+obj.x).^3);
        end
        
        function [u_out] = Apply_Misfit_Hessian(obj,u_in,u,z)
            u_out = obj.M*u_in;
        end
        
        function [u_opt] = Load_Optimal_u(obj)
            u_opt = (1+obj.x).^3;
        end
        
        function [z_opt] = Load_Optimal_z(obj)
            z_opt = 1+obj.x;
        end
        
        function [Z] = Load_Z_Data(obj)
            Z = zeros(obj.m,2);
            Z(:,1) = 1 + obj.x;
            Z(:,2) = obj.x + obj.x.^2;
        end
        
        function [D] = Load_d_Data(obj)
            Z = obj.Load_Z_Data();
            D = .2*(Z.^3);
        end
    end
    
    methods
        function obj = HDSA_MD_Interface_synthetic_test_with_hessian_gevp(m)
            obj.m = m;
            obj.x = linspace(0,1,m)';
            
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
            
            E_d = (1.e-6)*S + M;
            E_u = (2.0)*( (5.e-2)*S + M );
            E_z = (1.e2)*( (1.e-2)*S + M );
            
            obj.W_d = E_d'*linsolve(M,E_d);
            obj.W_u = E_u'*linsolve(M,E_u);
            obj.W_z = E_z'*linsolve(M,E_z);
        end
        
    end
    
end