classdef HDSA_MD_Update < handle
    
    properties
        md_interface;
        post_data;
        u_opt;
        z_opt;
        gevp;
    end
        
    methods
        function obj = HDSA_MD_Update(md_interface)
            obj.md_interface = md_interface;
            obj.post_data = HDSA_Bayes_Posterior_Data();
            obj.u_opt = obj.md_interface.Load_Optimal_u();
            obj.z_opt = obj.md_interface.Load_Optimal_z();
        end
        
        function [] = Compute_Posterior_Data(obj,alpha_d,num_samples)
            obj.post_data.Compute_Posterior_Data(obj.md_interface,alpha_d,obj.u_opt,obj.z_opt,num_samples);
        end
        
        function [] = Compute_Hessian_GEVP(obj,num_evals,oversampling)
            obj.gevp = Hessian_GEVP(obj.md_interface,obj.z_opt);
            obj.gevp.Compute_Hessian_GEVP(num_evals,oversampling);
        end
        
        function [delta_mean,delta_samples] = Posterior_Discrepancy_Samples(obj,z)
            m = size(obj.post_data.u_ell,1);
            p = size(z,2);
            delta_mean = cell(p,1);
            delta_samples = cell(p,1);
            
            Zc = obj.post_data.Z(:,2:end) - obj.z_opt;
            Zc_W_z_Inv_Zc = Zc'*(obj.post_data.W_z_inv_Z(:,2:end) - obj.post_data.W_z_inv_z_opt);
            
            for k = 1:p
                dz = z(:,k) - obj.z_opt;
                
                coeffs = 1 + (obj.post_data.W_z_inv_Z-obj.post_data.W_z_inv_z_opt)'*dz;
                delta_mean_k = obj.post_data.u_ell*coeffs;
                
                delta_samples_k = zeros(m,obj.post_data.num_samples);
                for i = 1:obj.post_data.N
                    sgi = sum(obj.post_data.g_vecs(:,i));
                    W_z_Inv_yi = obj.post_data.W_z_inv_Z*obj.post_data.g_vecs(:,i) - sgi*obj.post_data.W_z_inv_z_opt;
                    
                    coeffs = obj.post_data.b_i_ell(i,:)'*(sgi + W_z_Inv_yi'*dz);
                    delta_mean_k = delta_mean_k - obj.post_data.u_i_ell{i}*coeffs;
                    
                    coeff = (1/sqrt(obj.post_data.Mu(i,i)))*(sgi + W_z_Inv_yi'*dz);
                    delta_samples_k = delta_samples_k + coeff*obj.post_data.ui_hat{i};
                end
                delta_mean_k = (1/obj.post_data.alpha_d)*delta_mean_k;
                delta_samples_k = sqrt(obj.post_data.alpha_d)*delta_samples_k;
                
                W_z_Inv_dz = obj.md_interface.Apply_W_z_Inverse(dz);
                tmp = dz'*W_z_Inv_dz - W_z_Inv_dz'*Zc*linsolve(Zc_W_z_Inv_Zc,Zc'*W_z_Inv_dz);
                if tmp < -1.e-13
                   disp('Error in Posterior Discrepancy Samples') 
                end
                breve_coeff = sqrt(abs(tmp));
                delta_samples_k = delta_samples_k + breve_coeff*obj.post_data.u_breve;
                
                delta_mean{k} = delta_mean_k;
                delta_samples{k} = delta_samples_k + delta_mean_k;
                
            end
        end
            
        function [z_update_mean,z_update_samples] = Posterior_Update_Samples(obj)
            z_update_mean = obj.Posterior_Update_Mean();
            
            m = size(obj.post_data.ui_hat{1},1);
            n = length(obj.z_opt);
            
            Btheta_hat = zeros(n,obj.post_data.num_samples);
            u_tmp1 = zeros(m,obj.post_data.num_samples);
            for i = 1:obj.post_data.N
                sgi = sum(obj.post_data.g_vecs(:,i));
                coeff = sgi/sqrt(obj.post_data.Mu(i,i));
                u_tmp1 = u_tmp1 + coeff*obj.post_data.ui_hat{i};
                
                coeff = obj.post_data.state_grad'*obj.post_data.ui_hat{i};
                W_z_Inv_yi = obj.post_data.W_z_inv_Z*obj.post_data.g_vecs(:,i) - sgi*obj.post_data.W_z_inv_z_opt;  
                Btheta_hat = Btheta_hat + (1/sqrt(obj.post_data.Mu(i,i)))*W_z_Inv_yi*coeff;
            end
            tmp = obj.md_interface.Apply_Misfit_Hessian(u_tmp1,obj.u_opt,obj.z_opt);
            Btheta_hat = Btheta_hat + obj.md_interface.Apply_Solution_Operator_z_Jacobian_Transpose(tmp,obj.z_opt);
            Btheta_hat = sqrt(obj.post_data.alpha_d)*Btheta_hat;
            
            Zc = obj.post_data.Z(:,2:end) - obj.z_opt;
            W_z_Inv_Zc = obj.post_data.W_z_inv_Z(:,2:end) - obj.post_data.W_z_inv_z_opt;
            Zc_W_z_Inv_Zc = Zc'*W_z_Inv_Zc;
            tmp1 = Zc'*obj.post_data.zbreve;
            tmp2 = linsolve(Zc_W_z_Inv_Zc,tmp1);
            tmp3 = W_z_Inv_Zc*tmp2;
            coeff = sqrt(obj.post_data.state_grad_W_u_inv_state_grad);
            Btheta_breve = coeff*(obj.post_data.zbreve - tmp3);
            
            if ~isempty(obj.gevp)
                z_update_samples = z_update_mean -   obj.gevp.Apply_Projected_RS_Hessian_Inverse(Btheta_hat+Btheta_breve);
            else
                z_update_samples = z_update_mean - obj.md_interface.Apply_RS_Hessian_Inverse(Btheta_hat+Btheta_breve,obj.z_opt);
            end
        end
        
        function [z_update_mean] = Posterior_Update_Mean(obj)
            
            N = obj.post_data.N;
            u = 0*obj.u_opt;
            for ell = 1:N
                u = u + obj.post_data.u_ell(:,ell);
                for i = 1:N
                    u = u - obj.post_data.b_i_ell(i,ell)*sum(obj.post_data.g_vecs(:,i))*obj.post_data.u_i_ell{i}(:,ell);
                end
            end
            tmp1 = obj.md_interface.Apply_Misfit_Hessian(u,obj.u_opt,obj.z_opt);
            z_tmp = obj.md_interface.Apply_Solution_Operator_z_Jacobian_Transpose(tmp1,obj.z_opt);
            
            for ell = 1:N
                z_tmp = z_tmp + (obj.post_data.state_grad'*obj.post_data.u_ell(:,ell))*(obj.post_data.W_z_inv_Z(:,ell)-obj.post_data.W_z_inv_z_opt);
                for i = 1:N
                    z_tmp = z_tmp - obj.post_data.b_i_ell(i,ell)*(obj.post_data.state_grad'*obj.post_data.u_i_ell{i}(:,ell))*(obj.post_data.W_z_inv_Z*obj.post_data.g_vecs(:,i)-sum(obj.post_data.g_vecs(:,i))*obj.post_data.W_z_inv_z_opt);
                end
            end
            
            z_tmp = (1/obj.post_data.alpha_d)*z_tmp;
            if ~isempty(obj.gevp)
                z_pert = -obj.gevp.Apply_Projected_RS_Hessian_Inverse(z_tmp);
            else
                z_pert = -obj.md_interface.Apply_RS_Hessian_Inverse(z_tmp,obj.z_opt);
            end
            z_update_mean = obj.z_opt + z_pert;
        end
        
    end
    
end