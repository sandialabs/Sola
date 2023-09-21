classdef HDSA_MD_Interface < handle
    
    properties

    end
    
    methods (Abstract, Access = public)
        
        %% Pure virtual functions
        
        [u_out] = Apply_W_d(obj,u_in);
        
        [z_out] = Apply_W_z_Inverse(obj,z_in);
        
        [u_out] = Apply_W_u_Plus_scalar_W_d_Inverse(obj,u_in,beta);
        
        [u_out] = Apply_W_u_Inverse(obj,u_in);
        
        [z_out] = Apply_Solution_Operator_z_Jacobian_Transpose(obj,u_in,z);
        
        [z_out] = Apply_RS_Hessian(obj,z_in,z);
        
        [grad_u] = Misfit_Gradient(obj,u,z);
        
        [u_out] = Apply_Misfit_Hessian(obj,u_in,u,z);
        
        [u_opt] = Load_Optimal_u(obj);
         
        [z_opt] = Load_Optimal_z(obj);
        
        [Z] = Load_Z_Data(obj);
        
        [D] = Load_d_Data(obj);
        
    end
    
    methods
        function obj = HDSA_MD_Interface()

        end
        
        % Factorize W_u^{-1}=F*F^T, function gives u_out=F*u_in
        % This function must be implemented to enable posterior update sampling
        function [u_out] = Apply_W_u_Inverse_Factor(obj,u_in)
            % Default implementation
            u_out = 0;
            disp('Apply_W_u_Inverse_Factor must be implemented to use sampling algorithms')
        end
        
        % Factorize (W_u+scalar*W_d)^{-1}=F*F^T, function gives u_out=F*u_in
        % This function must be implemented to enable posterior update sampling
        function [u_out] = Apply_W_u_Plus_scalar_W_d_Inverse_Factor(obj,u_in,scalar)
            % Default implementation
            u_out = 0;
            disp('Apply_W_u_Plus_scalar_W_d_Inverse_Factor must be implemented to use sampling algorithms')
        end
        
        % Factorize W_z^{-1} = F*F^T, function gives z_out = F*z_in
        % This function must be implemented to enable posterior update sampling
        function [z_out] = Apply_W_z_Inverse_Factor(obj,z_in)
            % Default implementation
            z_out = 0;
            disp('Apply_W_z_Inverse_Factor must be implemented to use sampling algorithms')         
        end
        
        % Apply W_z matrix
        % This function must be implemented to enable Hessian GEVP
        function [z_out] = Apply_W_z(obj,z_in)
            % Default implementation
            z_out = 0;
            disp('Apply_W_z must be implemented to use Hessian GEVP')
        end
        
        function [z_out] = Apply_RS_Hessian_Inverse(obj,z_in,z)
            z_out = 0*z_in;
            for k = 1:size(z_in,2)
                tol = 1.e-7;
                max_iter = length(z);
                [z_out(:,k),flag,relres,iter,resvec] = pcg(@(x)obj.Apply_RS_Hessian(x,z),z_in(:,k),tol,max_iter);
                if flag~=0
                    disp('CG did not converge')
                end
            end
        end
        
    end
    
end