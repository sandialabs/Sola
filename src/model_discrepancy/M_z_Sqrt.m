classdef M_z_Sqrt < Matrix_Sqrt
    
    properties
        MD_Interface_Elliptic_Prior;
    end
    
    methods (Access = public)
        
        
        function [vec_out] = Matrix_Apply(this,vec_in)
           vec_out = this.MD_Interface_Elliptic_Prior.Apply_M_z(vec_in); 
        end
        
    end
    
    methods
        function this = M_z_Sqrt(MD_Interface_Elliptic_Prior)
            this.MD_Interface_Elliptic_Prior = MD_Interface_Elliptic_Prior;
        end
        
        % Overload this function with a preconditioner G such that G^T*G
        % approxiates the matrix implemented in the "Matrix_Apply" function
        function [vec_out] = Matrix_Preconditioner_Apply(this,vec_in)
            vec_out = vec_in;
        end
      
    end
    
end