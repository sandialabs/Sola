classdef Randomized_GEVP_Test < Randomized_GEVP
    
    properties
        A;
        M;
    end
    
    methods (Access = public)
        
        function [vec_out] = Apply_Operator(obj,vec_in)
            vec_out = linsolve(obj.A,vec_in);
        end
        
        function [vec_out] = Apply_Weighting_Operator(obj,vec_in)
            vec_out = obj.M*vec_in;
        end
        
        function [vec_out] = Apply_Weighting_Operator_Inverse(obj,vec_in)
            vec_out = linsolve(obj.M,vec_in);
        end
        
        function [vec_out] = Apply_Weighting_Operator_Inverse_Factor(obj,vec_in)
            vec_out = linsolve(sqrtm(obj.M),vec_in);
        end
        
    end
    
    methods
        function obj = Randomized_GEVP_Test(m)
            vec = zeros(m,1);
            obj@Randomized_GEVP(vec);
            
            h = 1/(m-1);
            
            M = diag(4*ones(1,m)) + diag(ones(1,m-1),1) + diag(ones(1,m-1),-1);
            M(1,1) = .5*M(1,1);
            M(end,end) = .5*M(end,end);
            M = (1/6)*h*M;
            obj.M = M;
            
            S = diag(2*ones(1,m)) + (-1)*diag(ones(1,m-1),1) + (-1)*diag(ones(1,m-1),-1);
            S(1,1) = .5*S(1,1);
            S(end,end) = .5*S(end,end);
            S = (1/h)*S;
            obj.A = (1.e-2)*S + M;
        end
           
    end
    
end