classdef HDSA_Abesa_MD_Interface_Elliptic_Prior_PDE_Test_Prob < HDSA_Abesa_MD_Interface_Elliptic_Prior

    properties
        E_u;
        E_z;
        E_d;
        M;
    end
    
    methods (Access = public)

        function [u_out] = Apply_E_u_Inverse(obj,u_in)
            u_out = linsolve(obj.E_u,u_in);
        end

        function [u_out] = Apply_E_u_Inverse_Transpose(obj,u_in)
            u_out = linsolve(obj.E_u',u_in);
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
      
        function [u_opt] = Load_Optimal_u(obj)
            u_opt = load('u_opt.mat').u_opt;
        end
        
        function [z_opt] = Load_Optimal_z(obj)
            z_opt = load('z_opt.mat').z_opt;
        end
        
        function [Z] = Load_Z_Data(obj)
            Z = load('Z.mat').Z;
        end
        
        function [D] = Load_d_Data(obj)
            D = load('D.mat').D;
        end
        
               
    end
    
    methods
        function obj = HDSA_Abesa_MD_Interface_Elliptic_Prior_PDE_Test_Prob(con_opt_obj,alpha_u,alpha_z)
            obj@HDSA_Abesa_MD_Interface_Elliptic_Prior(con_opt_obj,alpha_u,alpha_z);
            
            obj.E_d = (1.e-6)*con_opt_obj.S + con_opt_obj.M;
            obj.E_u = (1*10^-2)*con_opt_obj.S + con_opt_obj.M;
            obj.E_z = (10^-3)*con_opt_obj.S + con_opt_obj.M;
            obj.M = con_opt_obj.M;
            
            num_sing_vals = 200;
            oversampling = 0;
            num_subspace_iters = 1;
            u_vec = zeros(con_opt_obj.m,1);
            obj.Compute_Elliptic_GSVD(num_sing_vals,oversampling,num_subspace_iters,u_vec);
            
            
        end
        
    end
    
end