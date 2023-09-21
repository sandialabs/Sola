classdef HDSA_Abesa_MD_Interface_Elliptic_Prior < HDSA_MD_Interface_Elliptic_Prior
    
    properties
        con_opt_obj;
        z_current;
        u_current;
        hessian_data;
    end
    
    methods (Abstract, Access = public)
        
        %% Pure virtual functions
        
        [u_out] = Apply_E_u_Inverse(obj,u_in);

        [u_out] = Apply_E_u_Inverse_Transpose(obj,u_in);
        
        [u_out] = Apply_M_u(obj,u_in);
        
        [u_out] = Apply_M_u_Inverse(obj,u_in);
        
        [z_out] = Apply_E_z_Inverse(obj,z_in);
        
        [z_out] = Apply_E_z_Inverse_Transpose(obj,z_in);
        
        [z_out] = Apply_M_z(obj,z_in);
        
        [u_out] = Apply_E_d(obj,u_in);
        
        [u_out] = Apply_E_d_Transpose(obj,u_in);
                
        [u_opt] = Load_Optimal_u(obj);
        
        [z_opt] = Load_Optimal_z(obj);
        
        [Z] = Load_Z_Data(obj);
        
        [D] = Load_d_Data(obj);
        
               
    end
    
    methods
        function obj = HDSA_Abesa_MD_Interface_Elliptic_Prior(con_opt_obj,alpha_u,alpha_z)
            obj@HDSA_MD_Interface_Elliptic_Prior(alpha_u,alpha_z);
            obj.con_opt_obj = con_opt_obj;
            obj.z_current = obj.Load_Optimal_z();
            [~,~,obj.hessian_data] = obj.con_opt_obj.Jhat(obj.z_current);
            m = (length(obj.hessian_data)-length(obj.z_current))/2;
            obj.u_current = obj.hessian_data(1:m);
        end

        function [z_out] = Apply_Solution_Operator_z_Jacobian_Transpose(obj,u_in,z)
            if norm(z-obj.z_current)~=0
                [~,~,obj.hessian_data] = obj.con_opt_obj.Jhat(obj.z_current);
                obj.z_current = z;
                obj.u_current = obj.hessian_data(1:obj.con_opt_obj.m);
            end
            tmp = obj.con_opt_obj.c_u_Transpose_Inverse_Apply(u_in,obj.u_current,z);
            z_out = -obj.con_opt_obj.c_z_Transpose_Apply(tmp,obj.u_current,z);
        end
        
        function [z_out] = Apply_RS_Hessian(obj,z_in,z)
            if norm(z-obj.z_current)~=0
                [~,~,obj.hessian_data] = obj.con_opt_obj.Jhat(obj.z_current);
                obj.z_current = z;
                obj.u_current = obj.hessian_data(1:obj.con_opt_obj.m);
            end
            z_out = obj.con_opt_obj.Jhat_hessVec(obj.hessian_data,z_in);
        end
        
        function [grad_u] = Misfit_Gradient(obj,u,z)
           [~,grad_u,~] = obj.con_opt_obj.Objective(u,z);
        end
        
        function [u_out] = Apply_Misfit_Hessian(obj,u_in,u,z)
            u_out = obj.con_opt_obj.Objective_uu_Apply(u_in,u,z);
        end
        
    end
    
end