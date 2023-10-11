classdef HDSA_MD_Interface_Elliptic_Prior < HDSA_MD_Interface
    
    properties
        sing_vecs_input;
        sing_vecs_output;
        sing_vals;
        alpha_u;
        alpha_z;
    end
    
    methods (Abstract, Access = public)
        
        %% Pure virtual functions

        [u_out] = Apply_E_u_Inverse(this,u_in);

        [u_out] = Apply_E_u_Inverse_Transpose(this,u_in);
        
        [u_out] = Apply_M_u(this,u_in);
        
        [u_out] = Apply_M_u_Inverse(this,u_in);
        
        [z_out] = Apply_E_z_Inverse(this,z_in);
        
        [z_out] = Apply_E_z_Inverse_Transpose(this,z_in);
        
        [z_out] = Apply_M_z(this,z_in);
        
        [u_out] = Apply_E_d(this,u_in);
        
        [u_out] = Apply_E_d_Transpose(this,u_in);
        
        [z_out] = Apply_Solution_Operator_z_Jacobian_Transpose(this,u_in,z);
        
        [z_out] = Apply_RS_Hessian(this,z_in,z);
        
        [grad_u] = Misfit_Gradient(this,u,z);
        
        [u_out] = Apply_Misfit_Hessian(this,u_in,u,z);
        
        [u_opt] = Load_Optimal_u(this);
        
        [z_opt] = Load_Optimal_z(this);
        
        [Z] = Load_Z_Data(this);
        
        [D] = Load_d_Data(this);
           
    end
    
    methods
        function this = HDSA_MD_Interface_Elliptic_Prior(alpha_u,alpha_z)
            this.alpha_u = alpha_u;
            this.alpha_z = alpha_z;
        end

        % This function must be implemented to enable Hessian GEVP
        function [z_out] = Apply_E_z(this,z_in)
            % Default implementation
            z_out = 0;
            disp('Apply_E_z must be implemented to use Hessian GEVP')
        end
        
        % This function must be implemented to enable Hessian GEVP
        function [z_out] = Apply_E_z_Transpose(this,z_in)
            % Default implementation
            z_out = 0;
            disp('Apply_E_z_Transpose must be implemented to use Hessian GEVP')
        end
        
        % This function must be implemented to enable Hessian GEVP
        function [z_out] = Apply_M_z_Inverse(this,z_in)
            % Default implementation
            z_out = 0;
            disp('Apply_M_z_Inverse must be implemented to use Hessian GEVP')
        end
        
        function [u_out] = Apply_W_d(this,u_in)
            tmp1 = this.Apply_E_d(u_in);
            tmp2 = this.Apply_M_u_Inverse(tmp1);
            u_out = this.Apply_E_d_Transpose(tmp2);
        end
        
        function [z_out] = Apply_W_z_Inverse(this,z_in)
            tmp1 = this.Apply_E_z_Inverse_Transpose(z_in);
            tmp2 = this.Apply_M_z(tmp1);
            z_out = this.alpha_z*this.Apply_E_z_Inverse(tmp2);
        end
        
        % Apply W_z
        % This function must be implemented to enable Hessian GEVP
        function [z_out] = Apply_W_z(this,z_in)
            tmp1 = this.Apply_E_z(z_in);
            tmp2 = this.Apply_M_z_Inverse(tmp1);
            z_out = (1/this.alpha_z)*this.Apply_E_z_Transpose(tmp2);
        end
        
        % Factorize W_z^{-1} = F*F^T, function gives z_out = F*z_in
        % This function must be implemented to enable posterior update sampling
        function [z_out] = Apply_W_z_Inverse_Factor(this,z_in)
            n = size(z_in,1);
            d = size(z_in,2);
            tmp = zeros(n,d);
            A = @(v) this.Apply_M_z(v);
            G = speye(n);
            tol = 1.e-8;
            for k = 1:d
                tmp(:,k) = krylov_sqrt(A,G,z_in(:,k),n,tol);
            end
            z_out = sqrt(this.alpha_z)*this.Apply_E_z_Inverse(tmp);
        end
        
        function [u_out] = Apply_W_u_Plus_scalar_W_d_Inverse(this,u_in,scalar)
            K = (this.sing_vals.^2)./(1+this.alpha_u*scalar*this.sing_vals.^2);
            u_out = this.alpha_u*this.sing_vecs_output*diag(K)*this.sing_vecs_output'*u_in;
        end
        
        function [u_out] = Apply_W_u_Inverse(this,u_in)
           u_out = this.alpha_u*this.sing_vecs_output*diag(this.sing_vals.^2)*this.sing_vecs_output'*u_in;
        end
        
        % Factorize W_u^{-1}=F*F^T, function gives u_out=F*u_in
        % This function must be implemented to enable posterior update sampling
        function [u_out] = Apply_W_u_Inverse_Factor(this,u_in)
            r = length(this.sing_vals);
            u_out = sqrt(this.alpha_u)*this.sing_vecs_output*diag(this.sing_vals)*u_in(1:r,:);
        end
        
        % Factorize (W_u+scalar*W_d)^{-1}=F*F^T, function gives u_out=F*u_in
        % This function must be implemented to enable posterior update sampling
        function [u_out] = Apply_W_u_Plus_scalar_W_d_Inverse_Factor(this,u_in,scalar)
            K = (this.sing_vals.^2)./(1+this.alpha_u*scalar*this.sing_vals.^2);
            r = length(this.sing_vals);
            u_out = sqrt(this.alpha_u)*this.sing_vecs_output*diag(sqrt(K))*u_in(1:r,:);
        end
        
        function [] = Compute_Elliptic_GSVD(this, num_sing_vals, oversampling, num_subspace_iters, u_vec)
            elliptic_gsvd = Elliptic_GSVD(this,u_vec,u_vec);
            [this.sing_vecs_input,this.sing_vecs_output,this.sing_vals] = elliptic_gsvd.Compute_GSVD(num_sing_vals, oversampling, num_subspace_iters);
        end
        
    end
    
end