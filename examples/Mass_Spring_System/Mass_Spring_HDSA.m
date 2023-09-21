classdef Mass_Spring_HDSA < HDSA_Abesa_MD_Interface_Elliptic_Prior
    
    properties
        E_u;
        E_z;
        E_d;
        P_z;
        M;
        H;
        evecs;
        evals;
    end
    
    methods
        function obj = Mass_Spring_HDSA(con_opt_obj,alpha_u,alpha_z)
            obj@HDSA_Abesa_MD_Interface_Elliptic_Prior(con_opt_obj,alpha_u,alpha_z);
            
            N = con_opt_obj.N;
            h = con_opt_obj.t_mesh(2)-con_opt_obj.t_mesh(1);
            M = diag(4*ones(1,N)) + diag(ones(1,N-1),1) + diag(ones(1,N-1),-1);
            M(1,1) = .5*M(1,1);
            M(end,end) = .5*M(end,end);
            M = (1/6)*h*M;
            obj.M = M;
            
            S = diag(2*ones(1,N)) + (-1)*diag(ones(1,N-1),1) + (-1)*diag(ones(1,N-1),-1);
            S(1,1) = .5*S(1,1);
            S(end,end) = .5*S(end,end);
            S = (1/h)*S;
            
            obj.P_z = con_opt_obj.P_z;
            
            I = eye(2);
            I(2,2) = 4;
            obj.E_u = (5.e-2)*kron(S,I) + kron(obj.M,I);
            obj.E_u(1:2,:) = 0;
            obj.E_u(1,1) = 10;
            obj.E_u(2,2) = 40;
            
            obj.E_z = (1.e-1)*S + obj.M;
            
            obj.E_d = (1.e-8)*kron(S,I) + kron(obj.M,I);
            
            num_sing_vals = 100;
            oversampling = 0;
            num_subspace_iters = 1;
            u_vec = zeros(con_opt_obj.m*N,1);
            obj.Compute_Elliptic_GSVD(num_sing_vals,oversampling,num_subspace_iters,u_vec);
        end
        
        function [u_out] = Apply_E_u_Inverse(obj,u_in)
            u_out = linsolve(obj.E_u,u_in);
        end
        
        function [u_out] = Apply_E_u_Inverse_Transpose(obj,u_in)
            u_out = linsolve(obj.E_u',u_in);
        end
        
        function [u_out] = Apply_M_u(obj,u_in)
            u_out = kron(eye(obj.con_opt_obj.m),obj.M)*u_in;
        end
        
        function [u_out] = Apply_M_u_Inverse(obj,u_in)
            u_out = linsolve(kron(eye(obj.con_opt_obj.m),obj.M),u_in);
        end
        
        function [z_out] = Apply_E_z_Inverse(obj,z_in)
            z_out = linsolve(obj.P_z'*obj.E_z*obj.P_z,z_in);
        end
        
        function [z_out] = Apply_E_z_Inverse_Transpose(obj,z_in)
            z_out = linsolve(obj.P_z'*obj.E_z'*obj.P_z,z_in);
        end
        
        function [z_out] = Apply_M_z(obj,z_in)
            z_out = obj.P_z'*obj.M*obj.P_z*z_in;
        end
        
        % This function must be implemented to enable Hessian GEVP
        function [z_out] = Apply_E_z(obj,z_in)
            z_out = obj.P_z'*obj.E_z*obj.P_z*z_in;
        end
        
        % This function must be implemented to enable Hessian GEVP
        function [z_out] = Apply_E_z_Transpose(obj,z_in)
            z_out = obj.P_z'*obj.E_z'*obj.P_z*z_in;
        end
        
        % This function must be implemented to enable Hessian GEVP
        function [z_out] = Apply_M_z_Inverse(obj,z_in)
            z_out = linsolve(obj.P_z'*obj.M*obj.P_z,z_in);
        end
        
        function [u_out] = Apply_E_d(obj,u_in)
            u_out = obj.E_d*u_in;
        end
        
        function [u_out] = Apply_E_d_Transpose(obj,u_in)
            u_out = obj.E_d'*u_in;
        end
                
        function [u_opt] = Load_Optimal_u(obj)
            u_opt = load('Optimization_Results.mat').u_lofi;
        end
        
        function [z_opt] = Load_Optimal_z(obj)
            z_opt = load('Optimization_Results.mat').z_lofi;
        end
        
        function [Z] = Load_Z_Data(obj)
            Z = load('Optimization_Results.mat').Z;
        end
        
        function [D] = Load_d_Data(obj)
            D = load('Optimization_Results.mat').D;
        end
    end
    
end