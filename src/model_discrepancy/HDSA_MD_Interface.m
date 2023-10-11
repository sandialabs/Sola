classdef HDSA_MD_Interface < handle
    
    properties

    end
    
    methods (Abstract, Access = public)
        
        %% Pure virtual functions
        
        [u_out] = Apply_W_d(this,u_in);
        
        [z_out] = Apply_W_z_Inverse(this,z_in);
        
        [u_out] = Apply_W_u_Plus_scalar_W_d_Inverse(this,u_in,beta);
        
        [u_out] = Apply_W_u_Inverse(this,u_in);
        
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
        function this = HDSA_MD_Interface()

        end
        
        % Factorize W_u^{-1}=F*F^T, function gives u_out=F*u_in
        % This function must be implemented to enable posterior update sampling
        function [u_out] = Apply_W_u_Inverse_Factor(this,u_in)
            % Default implementation
            u_out = 0;
            disp('Apply_W_u_Inverse_Factor must be implemented to use sampling algorithms')
        end
        
        % Factorize (W_u+scalar*W_d)^{-1}=F*F^T, function gives u_out=F*u_in
        % This function must be implemented to enable posterior update sampling
        function [u_out] = Apply_W_u_Plus_scalar_W_d_Inverse_Factor(this,u_in,scalar)
            % Default implementation
            u_out = 0;
            disp('Apply_W_u_Plus_scalar_W_d_Inverse_Factor must be implemented to use sampling algorithms')
        end
        
        % Factorize W_z^{-1} = F*F^T, function gives z_out = F*z_in
        % This function must be implemented to enable posterior update sampling
        function [z_out] = Apply_W_z_Inverse_Factor(this,z_in)
            % Default implementation
            z_out = 0;
            disp('Apply_W_z_Inverse_Factor must be implemented to use sampling algorithms')         
        end
        
        % Apply W_z matrix
        % This function must be implemented to enable Hessian GEVP
        function [z_out] = Apply_W_z(this,z_in)
            % Default implementation
            z_out = 0;
            disp('Apply_W_z must be implemented to use Hessian GEVP')
        end
        
        function [z_out] = Apply_RS_Hessian_Inverse(this,z_in,z)
            z_out = 0*z_in;
            for k = 1:size(z_in,2)
                tol = 1.e-7;
                max_iter = length(z);
                [z_out(:,k),flag,relres,iter,resvec] = pcg(@(x)this.Apply_RS_Hessian(x,z),z_in(:,k),tol,max_iter);
                if flag~=0
                    disp('CG did not converge')
                end
            end
        end
        
    end
    
end