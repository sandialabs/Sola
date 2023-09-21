classdef HDSA_Bayes_Posterior_Data < handle
    
    properties
        alpha_d;
        N;
        Z;
        D;
        W_z_inv_Z;
        W_z_inv_z_opt;
        G;
        g_vecs;
        Mu;
        state_grad;
        u_ell;
        u_i_ell;
        a_ell;
        b_i_ell;
        
        num_samples;
        ui_hat;
        u_breve;
        state_grad_W_u_inv_state_grad;
        zbreve; 
    end
    
    methods
        function obj = HDSA_Bayes_Posterior_Data()

        end
        
        function [] = Compute_Posterior_Data(obj, md_interface, alpha_d_in, u_opt, z_opt, num_samples)
            obj.alpha_d = alpha_d_in;
            obj.num_samples = num_samples;
            obj.Z = md_interface.Load_Z_Data();
            obj.D = md_interface.Load_d_Data();
            obj.N = size(obj.D,2);
            obj.state_grad = md_interface.Misfit_Gradient(u_opt,z_opt);
            
            obj.W_z_inv_Z = md_interface.Apply_W_z_Inverse(obj.Z);
            obj.W_z_inv_z_opt = md_interface.Apply_W_z_Inverse(z_opt);
            obj.G = (1+obj.W_z_inv_z_opt'*z_opt) - obj.Z'*obj.W_z_inv_z_opt - obj.W_z_inv_z_opt'*obj.Z + obj.Z'*obj.W_z_inv_Z;
            [obj.g_vecs,obj.Mu] = eig(obj.G);
            
            W_d_Y = md_interface.Apply_W_d(obj.D);
            obj.u_ell = md_interface.Apply_W_u_Inverse(W_d_Y);
            obj.u_i_ell = cell(obj.N,1);
            W_d_u_ell = md_interface.Apply_W_d(obj.u_ell);
            for i = 1:obj.N
                obj.u_i_ell{i} = (1/obj.alpha_d)*md_interface.Apply_W_u_Plus_scalar_W_d_Inverse(W_d_u_ell,obj.Mu(i,i)/obj.alpha_d);
            end
            
            obj.a_ell = zeros(obj.N,1);
            obj.b_i_ell = zeros(obj.N,obj.N);
            for ell = 1:obj.N
                obj.a_ell(ell) = 1 - obj.W_z_inv_z_opt'*(obj.Z(:,ell)-z_opt);
                for i = 1:obj.N
                    obj.b_i_ell(i,ell) = (obj.Z*obj.g_vecs(:,i))'*(obj.W_z_inv_Z(:,ell)-obj.W_z_inv_z_opt) + sum(obj.g_vecs(:,i))*obj.a_ell(ell);
                end
            end
            
            if obj.num_samples > 0
                
                obj.ui_hat = cell(obj.N,1);
                m = size(obj.u_ell,1);
                for i = 1:obj.N
                    Omega = randn(m,obj.num_samples);
                    obj.ui_hat{i} = (1/sqrt(obj.alpha_d))*md_interface.Apply_W_u_Plus_scalar_W_d_Inverse_Factor(Omega,obj.Mu(i,i)/obj.alpha_d);
                end
            
                Omega = randn(m,obj.num_samples);
                obj.u_breve = md_interface.Apply_W_u_Inverse_Factor(Omega);
                
                obj.state_grad_W_u_inv_state_grad = md_interface.Apply_W_u_Inverse(obj.state_grad)'*obj.state_grad;
                n = length(z_opt);
                Omega = randn(n,obj.num_samples);
                obj.zbreve = md_interface.Apply_W_z_Inverse_Factor(Omega);
                
            end
            
        end
    
        
    end
    
end