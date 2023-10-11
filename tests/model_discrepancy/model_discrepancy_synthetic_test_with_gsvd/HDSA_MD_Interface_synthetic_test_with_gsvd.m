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
        
        function [u_out] = Apply_E_u_Inverse(this,u_in)
            u_out = linsolve(this.E_u,u_in);
        end
        
        function [u_out] = Apply_E_u_Inverse_Transpose(this,u_in)
            u_out = linsolve(this.E_u,u_in);
        end
        
        function [u_out] = Apply_M_u(this,u_in)
            u_out = this.M*u_in;
        end
        
        function [u_out] = Apply_M_u_Inverse(this,u_in)
            u_out = linsolve(this.M,u_in);
        end

        function [z_out] = Apply_E_z_Inverse(this,z_in)
            z_out = linsolve(this.E_z,z_in);
        end
        
        function [z_out] = Apply_E_z_Inverse_Transpose(this,z_in)
            z_out = linsolve(this.E_z',z_in);
        end
        
        function [z_out] = Apply_M_z(this,z_in)
            z_out = this.M*z_in;
        end
        
        function [u_out] = Apply_E_d(this,u_in)
            u_out = this.E_d*u_in;
        end
        
        function [u_out] = Apply_E_d_Transpose(this,u_in)
            u_out = this.E_d'*u_in;
        end
                
        function [z_out] = Apply_Solution_Operator_z_Jacobian_Transpose(this,u_in,z)
            z_out = 3*diag(z.^2)*u_in;
        end
            
        % This implementation assumes that it is evaluated at the optimal z so that 
        % the adjoint=0, a more general impl a term multiplied by the adjoint variable
        function [z_out] = Apply_RS_Hessian(this,z_in,z)
            tmp1 = 3*diag(z.^2)*z_in;
            tmp2 = this.M*tmp1;
            z_out = 3*diag(z.^2)*tmp2;
        end
        
        function [grad_u] = Misfit_Gradient(this,u,z)
            grad_u = this.M*(u-(1+this.x).^3);
        end
        
        function [u_out] = Apply_Misfit_Hessian(this,u_in,u,z)
            u_out = this.M*u_in;
        end
        
        function [u_opt] = Load_Optimal_u(this)
            u_opt = (1+this.x).^3;
        end
        
        function [z_opt] = Load_Optimal_z(this)
            z_opt = 1+this.x;
        end
        
        function [Z] = Load_Z_Data(this)
            Z = zeros(this.m,2);
            Z(:,1) = 1 + this.x;
            Z(:,2) = this.x + this.x.^2;
        end
        
        function [D] = Load_d_Data(this)
            Z = this.Load_Z_Data();
            D = .2*(Z.^3);
        end
    end
    
    methods
        function this = HDSA_MD_Interface_synthetic_test_with_gsvd(m,alpha_u,alpha_w)
            this@HDSA_MD_Interface_Elliptic_Prior(alpha_u,alpha_w);
            
            this.m = m;
            this.x = linspace(0,1,m)';
            
            h = this.x(2)-this.x(1);
            M = diag(4*ones(1,m)) + diag(ones(1,m-1),1) + diag(ones(1,m-1),-1);
            M(1,1) = .5*M(1,1);
            M(end,end) = .5*M(end,end);
            M = (1/6)*h*M;
            this.M = M;
            
            S = diag(2*ones(1,m)) + (-1)*diag(ones(1,m-1),1) + (-1)*diag(ones(1,m-1),-1);
            S(1,1) = .5*S(1,1);
            S(end,end) = .5*S(end,end);
            S = (1/h)*S;
            this.S = S;
            
            this.E_d = (1.e-6)*S + M;
            this.E_u = (5.e-2)*S + M;
            this.E_z = (1.e-2)*S + M;
            
            num_sing_vals = 50;
            oversampling = 1;
            num_subspace_iters = 2;
            u_vec = zeros(m,1);
            this.Compute_Elliptic_GSVD(num_sing_vals,oversampling,num_subspace_iters,u_vec);
            
        end
        
    end
    
end