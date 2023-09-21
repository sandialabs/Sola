classdef HDSA_MD_Interface_synthetic_test_with_gsvd < HDSA_MD_Interface_Elliptic_Prior
    
    properties
        m; % Mesh resolution
        x; % Mesh nodes on [0,1]
        S; % Stiffness matrix
        M; % Mass matrix
        E_d; % Discrepancy precision matrix elliptic operator
        E_u; % State precision matrix elliptic operator
        E_z; % Opt precision matrix elliptic operator
    end
    
    methods (Access = public)
        
        function [u_out] = Apply_E_u_Inverse(obj,u_in)
            u_out = linsolve(obj.E_u,u_in);
        end
        
        function [u_out] = Apply_E_u_Inverse_Transpose(obj,u_in)
            u_out = linsolve(obj.E_u,u_in);
        end
        
        function [u_out] = Apply_M_u(obj,u_in)
            u_out = obj.M*u_in;
        end
        
        function [u_out] = Apply_M_u_Inverse(obj,u_in)
            u_out = linsolve(obj.M,u_in);
        end

        function [z_out] = Apply_E_z_Inverse(obj,z_in)
            z_out = linsolve(obj.E_z,z_in);
        end
        
        function [z_out] = Apply_E_z_Inverse_Transpose(obj,z_in)
            z_out = linsolve(obj.E_z',z_in);
        end
        
        function [z_out] = Apply_M_z(obj,z_in)
            z_out = obj.M*z_in;
        end
        
        function [u_out] = Apply_E_d(obj,u_in)
            u_out = obj.E_d*u_in;
        end
        
        function [u_out] = Apply_E_d_Transpose(obj,u_in)
            u_out = obj.E_d'*u_in;
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
        function obj = HDSA_MD_Interface_synthetic_test_with_gsvd(m,alpha_u,alpha_w)
            obj@HDSA_MD_Interface_Elliptic_Prior(alpha_u,alpha_w);
            
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
            
            obj.E_d = (1.e-6)*S + M;
            obj.E_u = (5.e-2)*S + M;
            obj.E_z = (1.e-2)*S + M;
            
            num_sing_vals = 50;
            oversampling = 1;
            num_subspace_iters = 2;
            u_vec = zeros(m,1);
            obj.Compute_Elliptic_GSVD(num_sing_vals,oversampling,num_subspace_iters,u_vec);
            
        end
        
    end
    
end