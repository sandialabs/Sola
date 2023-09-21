classdef Randomized_GSVD < handle
    
    properties
        vec_in;
        vec_out;
    end
    
    methods (Abstract, Access = public)
        
        %% Pure virtual functions

        [vec_out] = Apply_Operator(obj,vec_in);

        [vec_out] = Apply_Operator_Transpose(obj,vec_in);
        
        [vec_out] = Apply_Input_Weighting_Operator(obj,vec_in);
        
        [vec_out] = Apply_Input_Weighting_Operator_Inverse(obj,vec_in);
        
        [vec_out] = Apply_Output_Weighting_Operator(obj,vec_in);
               
    end
    
    methods
        function obj = Randomized_GSVD(vec_in,vec_out)
            obj.vec_in = 0*vec_in;
            obj.vec_out = 0*vec_out;
        end
   
        function [sing_vecs_input,sing_vecs_output,sing_vals] = Compute_GSVD(obj, num_sing_vals, oversampling, num_subspace_iters)

            kpp = num_sing_vals + oversampling;
            m = length(obj.vec_in);
            
            Omega = randn(m,kpp);
            Y = obj.Apply_Operator(Omega);
            
            [Q,WQ] = obj.CholQR(Y,'output_weighting');
            
            for j = 1:num_subspace_iters
                Y_subspace_iter = obj.Apply_Operator_Transpose(WQ);
                [~,WQ_subspace_iter] = obj.CholQR(Y_subspace_iter,'input_weighting_inverse');
                Y = obj.Apply_Operator(WQ_subspace_iter);
                [Q,WQ] = obj.CholQR(Y,'output_weighting');
            end
            
            B = obj.Apply_Operator_Transpose(WQ);
            Tinv_B = obj.Apply_Input_Weighting_Operator_Inverse(B);
            
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
        
        function [Q,WQ] = CholQR(obj,Z,type)
            [Q_Z,~] = qr(Z,0);
            if strcmp(type,'input_weighting')
                W_Q_Z = obj.Apply_Input_Weighting_Operator(Q_Z);
            elseif strcmp(type,'input_weighting_inverse')
                W_Q_Z = obj.Apply_Input_Weighting_Operator_Inverse(Q_Z);
            elseif strcmp(type,'output_weighting')
                W_Q_Z = obj.Apply_Output_Weighting_Operator(Q_Z);
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