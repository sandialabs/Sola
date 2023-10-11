classdef Randomized_GEVP_Test < Randomized_GEVP
    
    properties
        A;
        M;
    end
    
    methods (Access = public)
        
        function [vec_out] = Apply_Operator(this,vec_in)
            vec_out = linsolve(this.A,vec_in);
        end
        
        function [vec_out] = Apply_Weighting_Operator(this,vec_in)
            vec_out = this.M*vec_in;
        end
        
        function [vec_out] = Apply_Weighting_Operator_Inverse(this,vec_in)
            vec_out = linsolve(this.M,vec_in);
        end
        
        function [vec_out] = Apply_Weighting_Operator_Inverse_Factor(this,vec_in)
            vec_out = linsolve(sqrtm(this.M),vec_in);
        end
        
    end
    
    methods
        function this = Randomized_GEVP_Test(m)
            vec = zeros(m,1);
            this@Randomized_GEVP(vec);
            
            h = 1/(m-1);
            
            M = diag(4*ones(1,m)) + diag(ones(1,m-1),1) + diag(ones(1,m-1),-1);
            M(1,1) = .5*M(1,1);
            M(end,end) = .5*M(end,end);
            M = (1/6)*h*M;
            this.M = M;
            
            S = diag(2*ones(1,m)) + (-1)*diag(ones(1,m-1),1) + (-1)*diag(ones(1,m-1),-1);
            S(1,1) = .5*S(1,1);
            S(end,end) = .5*S(end,end);
            S = (1/h)*S;
            this.A = (1.e-2)*S + M;
        end
           
    end
    
end