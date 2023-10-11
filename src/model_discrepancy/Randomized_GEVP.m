classdef Randomized_GEVP < handle
    
    properties
        vec;
    end
    
    methods (Abstract, Access = public)
        
        %% Pure virtual functions

        [vec_out] = Apply_Operator(this,vec_in);
        
        [vec_out] = Apply_Weighting_Operator(this,vec_in);
        
        [vec_out] = Apply_Weighting_Operator_Inverse(this,vec_in);
        
        [vec_out] = Apply_Weighting_Operator_Inverse_Factor(this,vec_in);
               
    end
    
    methods
        function this = Randomized_GEVP(vec)
            this.vec = 0*vec;
        end
           
        function [evecs,evals] = Compute_GEVP(this, num_evals, oversampling)

            kpp = num_evals + oversampling;
            m = length(this.vec);
            
            Omega = this.Apply_Weighting_Operator_Inverse_Factor(randn(m,kpp)); 
            tmp = this.Apply_Operator(Omega);
            Y = this.Apply_Weighting_Operator_Inverse(tmp);
            
            Q = this.CholQR(Y,'weighting');
            
            AQ = this.Apply_Operator(Q);
            T = Q'*AQ;
            
            R_T = chol(T);
            M = AQ*linsolve(R_T,eye(size(R_T,1)));
            [~,WQ,R] = this.CholQR(M,'weighting_inverse');
            [U_M,Sigma_M,~] = svd(R);
            
            evecs = WQ*U_M(:,1:num_evals);
            evals = diag(Sigma_M(1:num_evals,1:num_evals)).^2;
        end
        
        function [Q,WQ,R] = CholQR(this,Z,type)
            [Q_Z,R_Z] = qr(Z,0);
            if strcmp(type,'weighting')
                W_Q_Z = this.Apply_Weighting_Operator(Q_Z);
            elseif strcmp(type,'weighting_inverse')
                W_Q_Z = this.Apply_Weighting_Operator_Inverse(Q_Z);
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