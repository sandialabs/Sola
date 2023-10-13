classdef Randomized_GSVD < handle
    
    properties
        vec_in;
        vec_out;
    end
    
    methods (Abstract, Access = public)
        
        %% Pure virtual functions

        [vec_out] = Apply_Operator(this,vec_in);

        [vec_out] = Apply_Operator_Transpose(this,vec_in);
        
        [vec_out] = Apply_Input_Weighting_Operator_Inverse(this,vec_in);
        
        [vec_out] = Apply_Output_Weighting_Operator(this,vec_in);
               
    end
    
    methods
        function this = Randomized_GSVD(vec_in,vec_out)
            this.vec_in = 0*vec_in;
            this.vec_out = 0*vec_out;
        end
   
        function [sing_vecs_input,sing_vecs_output,sing_vals] = Compute_GSVD(this, num_sing_vals, oversampling, num_subspace_iters)

            kpp = num_sing_vals + oversampling;
            m = length(this.vec_in);
            
            Omega = randn(m,kpp);
            Y = this.Apply_Operator(Omega);
            
            [Q,WQ] = this.CholQR(Y,'output_weighting');
            
            for j = 1:num_subspace_iters
                Y_subspace_iter = this.Apply_Operator_Transpose(WQ);
                [~,WQ_subspace_iter] = this.CholQR(Y_subspace_iter,'input_weighting_inverse');
                Y = this.Apply_Operator(WQ_subspace_iter);
                [Q,WQ] = this.CholQR(Y,'output_weighting');
            end
            
            B = this.Apply_Operator_Transpose(WQ);
            Tinv_B = this.Apply_Input_Weighting_Operator_Inverse(B);
            
            C = B'*Tinv_B;
            R_B = chol(C);
            R_B_inv = linsolve(R_B,eye(size(R_B,1)));
            Q_B = Tinv_B*R_B_inv;
            
            [U,Sigma,V] = svd(R_B');
            
            sing_vecs_input = Q_B*V;
            sing_vecs_output = Q*U;
            sing_vals = diag(Sigma);
            
            sing_vecs_input = sing_vecs_input(:,1:num_sing_vals);
            sing_vecs_output = sing_vecs_output(:,1:num_sing_vals);
            sing_vals = sing_vals(1:num_sing_vals);
        end
        
        function [Q,WQ] = CholQR(this,Z,type)
            [Q_Z,~] = qr(Z,0);
            if strcmp(type,'input_weighting_inverse')
                W_Q_Z = this.Apply_Input_Weighting_Operator_Inverse(Q_Z);
            elseif strcmp(type,'output_weighting')
                W_Q_Z = this.Apply_Output_Weighting_Operator(Q_Z);
            else
                disp('Error specifying type in CholQR')
            end
            
            C = Q_Z'*W_Q_Z;
            R_C = chol(C);
            
            R_C_inv = linsolve(R_C,eye(size(R_C,1)));
            Q = Q_Z*R_C_inv;
            WQ = W_Q_Z*R_C_inv;
        end
        
    end
    
end