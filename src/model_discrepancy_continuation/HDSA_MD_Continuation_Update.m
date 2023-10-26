classdef HDSA_MD_Continuation_Update < handle
    
    properties
        md_continuation_interface;
        post_data;
        u_opt;
        z_opt;
        num_continuation_steps;
        step_size;
        
        W_z_inv_Z_minus_z_opt;
        W_z_inv_yi;
        si;
    end
    
    methods
        function this = HDSA_MD_Continuation_Update(md_continuation_interface,num_continuation_steps)
            this.md_continuation_interface = md_continuation_interface;
            this.post_data = HDSA_Bayes_Posterior_Data();
            this.u_opt = this.md_continuation_interface.md_interface.Load_Optimal_u();
            this.z_opt = this.md_continuation_interface.md_interface.Load_Optimal_z();
            this.num_continuation_steps = num_continuation_steps;
            this.step_size = 1/num_continuation_steps;
        end
        
        function [] = Compute_Posterior_Data(this,alpha_d,num_samples)
            this.post_data.Compute_Posterior_Data(this.md_continuation_interface.md_interface,alpha_d,this.u_opt,this.z_opt,num_samples);
            this.W_z_inv_Z_minus_z_opt = this.post_data.W_z_inv_Z - this.post_data.W_z_inv_z_opt;
            this.W_z_inv_yi = 0*this.W_z_inv_Z_minus_z_opt;
            for i = 1:this.post_data.N
                this.W_z_inv_yi(:,i) = this.post_data.W_z_inv_Z*this.post_data.g_vecs(:,i) - sum(this.post_data.g_vecs(:,i))*this.post_data.W_z_inv_z_opt;
                this.si(i) = sum(this.post_data.g_vecs(:,i)) - this.z_opt'*this.W_z_inv_yi(:,i);
            end
        end
                
        function [u,z] = Posterior_Update_Mean(this)
            u = zeros(length(this.u_opt),this.num_continuation_steps+1);
            z = zeros(length(this.z_opt),this.num_continuation_steps+1);
            t = linspace(0,1,this.num_continuation_steps+1);
            
            u(:,1) = this.u_opt;
            z(:,1) = this.z_opt;
            
            for k = 1:this.num_continuation_steps
                Btheta_n = this.Apply_B(u(:,k),z(:,k),t(k));
                z_pert = -this.Apply_Parameterized_RS_Hessian_Inverse(Btheta_n,u(:,k),z(:,k),t(k));
                z(:,k+1) = z(:,k) + this.step_size*z_pert;
                u(:,k+1) = this.md_continuation_interface.State_Solve(z(:,k+1));
            end
        end
        
        function [Btheta_n] = Apply_B(this,u_n,z_n,t_n)
            u_tmp1 = this.Apply_Discrepancy_theta_Jacobian(z_n);
            u_tmp2 = this.md_continuation_interface.md_interface.Apply_Misfit_Hessian(u_tmp1,u_n,z_n);
            z_tmp1 = this.md_continuation_interface.md_interface.Apply_Solution_Operator_z_Jacobian_Transpose(u_tmp2,z_n);
            
            z_tmp2 = this.Apply_Discrepancy_z_Jacobian_transpose(u_tmp2,t_n);
            
            state_grad = this.md_continuation_interface.md_interface.Misfit_Gradient(u_n,z_n);
            z_tmp3 = this.Apply_Discrepancy_z_theta_Hessian(state_grad);
            
            Btheta_n = z_tmp1 + z_tmp2 + z_tmp3;
        end

        function [z_out] = Apply_Parameterized_RS_Hessian_Inverse(this,z_in,u_n,z_n,t_n)
            z_out = 0*z_in;
            for k = 1:size(z_in,2)
                tol = 1.e-7;
                max_iter = length(z_n);
                [z_out(:,k),flag,relres,iter,resvec] = pcg(@(x)this.Apply_Parameterized_RS_Hessian(x,u_n,z_n,t_n),z_in(:,k),tol,max_iter);
                if flag~=0
                    disp('CG did not converge')
                end
            end
        end
        
        function [z_out] = Apply_Parameterized_RS_Hessian(this,z_in,u_n,z_n,t_n)
            z_out = this.md_continuation_interface.md_interface.Apply_RS_Hessian(z_in,z_n);
            
            u_tmp1 = this.Apply_Discrepancy_z_Jacobian(z_in,t_n);
            u_tmp2 = this.md_continuation_interface.md_interface.Apply_Misfit_Hessian(u_tmp1,u_n,z_n);
            z_out = z_out + this.md_continuation_interface.md_interface.Apply_Solution_Operator_z_Jacobian_Transpose(u_tmp2,z_n);
            
            z_out = z_out + this.Apply_Discrepancy_z_Jacobian_transpose(u_tmp2,t_n);
            
            u_tmp3 = this.md_continuation_interface.Apply_Solution_Operator_Jacobian(z_in,z_n);
            u_tmp4 = this.md_continuation_interface.md_interface.Apply_Misfit_Hessian(u_tmp3,u_n,z_n);
            z_out = z_out + this.Apply_Discrepancy_z_Jacobian_transpose(u_tmp4,t_n);
        end
        
        function [u_out] = Apply_Discrepancy_z_Jacobian(this,z_in,t_n)
            N = this.post_data.N;
            u = 0*this.u_opt;
            for ell = 1:N
                u = u + ((this.post_data.W_z_inv_Z(:,ell)-this.post_data.W_z_inv_z_opt)'*z_in)*this.post_data.u_ell(:,ell);
                for i = 1:N
                    u = u - this.post_data.b_i_ell(i,ell)*((this.post_data.W_z_inv_Z*this.post_data.g_vecs(:,i)-sum(this.post_data.g_vecs(:,i))*this.post_data.W_z_inv_z_opt)'*z_in)*this.post_data.u_i_ell{i}(:,ell);
                end
            end

            u_out = t_n*(1/this.post_data.alpha_d)*u;
        end
        
        function [z_out] = Apply_Discrepancy_z_Jacobian_transpose(this,u_in,t_n) 
            N = this.post_data.N;
            z = 0*this.z_opt;
            for ell = 1:N
                z = z + (this.post_data.u_ell(:,ell)'*u_in)*(this.post_data.W_z_inv_Z(:,ell)-this.post_data.W_z_inv_z_opt);
                for i = 1:N
                    z = z - this.post_data.b_i_ell(i,ell)*(this.post_data.u_i_ell{i}(:,ell)'*u_in)*(this.post_data.W_z_inv_Z*this.post_data.g_vecs(:,i)-sum(this.post_data.g_vecs(:,i))*this.post_data.W_z_inv_z_opt);
                end
            end

            z_out = t_n*(1/this.post_data.alpha_d)*z;
        end
        
        function [u_out] = Apply_Discrepancy_theta_Jacobian(this,z_n)
            N = this.post_data.N;
            u_out = 0*this.u_opt;
            for ell = 1:N
                coeff = this.post_data.a_ell(ell) + z_n'*this.W_z_inv_Z_minus_z_opt(:,ell);
                u_out = u_out + coeff*this.post_data.u_ell(:,ell);
                for i = 1:N
                    coeff = this.post_data.b_i_ell(i,ell)*(this.si(i) + z_n'*this.W_z_inv_yi(:,i));
                    u_out = u_out - coeff*this.post_data.u_i_ell{i}(:,ell);
                end
            end
            u_out = (1/this.post_data.alpha_d)*u_out;
        end
        
        function [z_out] = Apply_Discrepancy_z_theta_Hessian(this,u)
            N = this.post_data.N;
            z_out = 0*this.z_opt;
            for ell = 1:N
                z_out = z_out + (u'*this.post_data.u_ell(:,ell))*(this.post_data.W_z_inv_Z(:,ell)-this.post_data.W_z_inv_z_opt);
                for i = 1:N
                    z_out = z_out - this.post_data.b_i_ell(i,ell)*(u'*this.post_data.u_i_ell{i}(:,ell))*(this.post_data.W_z_inv_Z*this.post_data.g_vecs(:,i)-sum(this.post_data.g_vecs(:,i))*this.post_data.W_z_inv_z_opt);
                end
            end
        end
        
    end
    
end