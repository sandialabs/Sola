classdef Mass_Matrix_Sqrt < Matrix_Sqrt
    
    properties
        inf_dim_prior;
    end
    
    methods (Access = public)
        
        
        function [vec_out] = Matrix_Apply(this,vec_in)
           vec_out = this.inf_dim_prior.Mass_Matrix_Apply(vec_in); 
        end
        
    end
    
    methods
        function this = Mass_Matrix_Sqrt(inf_dim_prior)
            this.inf_dim_prior = inf_dim_prior;
        end
        
      
    end
    
end