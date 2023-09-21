classdef Randomized_GEVP < handle
    
    properties
        vec;
    end
    
    methods (Abstract, Access = public)
        
        %% Pure virtual functions

        [vec_out] = Apply_Operator(obj,vec_in);
        
        [vec_out] = Apply_Weighting_Operator(obj,vec_in);
        
        [vec_out] = Apply_Weighting_Operator_Inverse(obj,vec_in);
        
        [vec_out] = Apply_Weighting_Operator_Inverse_Factor(obj,vec_in);
               
    end
    
    methods
        function obj = Randomized_GEVP(vec)
            obj.vec = 0*vec;
        end
           
        function [evecs,evals] = Compute_GEVP(obj, num_evals, oversampling)

            kpp = num_evals + oversampling;
            m = length(obj.vec);
            
            Omega = obj.Apply_Weighting_Operator_Inverse_Factor(randn(m,kpp)); 
            tmp = obj.Apply_Operator(Omega);
            Y = obj.Apply_Weighting_Operator_Inverse(tmp);
            
            Q = obj.CholQR(Y,'weighting');
            
            AQ = obj.Apply_Operator(Q);
            T = Q'*AQ;
            
            R_T = chol(T);
            M = AQ*linsolve(R_T,eye(size(R_T,1)));
            [~,WQ,R] = obj.CholQR(M,'weighting_inverse');
            [U_M,Sigma_M,~] = svd(R);
            
            evecs = WQ*U_M(:,1:num_evals);
            evals = diag(Sigma_M(1:num_evals,1:num_evals)).^2;
        end
        
        function [Q,WQ,R] = CholQR(obj,Z,type)
            [Q_Z,R_Z] = qr(Z,0);
            if strcmp(type,'weighting')
                W_Q_Z = obj.Apply_Weighting_Operator(Q_Z);
            elseif strcmp(type,'weighting_inverse')
                W_Q_Z = obj.Apply_Weighting_Operator_Inverse(Q_Z);
            else
                disp('Error specifying type in CholQR')
            end
            C = Q_Z'*W_Q_Z;
            R_C = chol(C);
            R_C_inv = linsolve(R_C,eye(size(R_C,1)));
            Q = Q_Z*R_C_inv;
            WQ = W_Q_Z*R_C_inv;
            R = R_C*R_Z;
        end
        
    end
    
end