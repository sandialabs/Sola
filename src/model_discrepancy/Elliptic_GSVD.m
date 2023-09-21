classdef Elliptic_GSVD < Randomized_GSVD
    
    properties
        md_interface_elliptic_prior;
    end
    
    methods (Access = public)
        
        function [vec_out] = Apply_Operator(obj,vec_in)
            vec_out = obj.md_interface_elliptic_prior.Apply_E_u_Inverse(vec_in);
        end

        function [vec_out] = Apply_Operator_Transpose(obj,vec_in)
            vec_out = obj.md_interface_elliptic_prior.Apply_E_u_Inverse_Transpose(vec_in);
        end
        
        function [vec_out] = Apply_Input_Weighting_Operator(obj,vec_in)
            vec_out = obj.md_interface_elliptic_prior.Apply_M_u_Inverse(vec_in);
        end
        
        function [vec_out] = Apply_Input_Weighting_Operator_Inverse(obj,vec_in)
            vec_out = obj.md_interface_elliptic_prior.Apply_M_u(vec_in);
        end
        
        function [vec_out] = Apply_Output_Weighting_Operator(obj,vec_in)
            vec_out = obj.md_interface_elliptic_prior.Apply_W_d(vec_in);
        end
               
    end
    
    methods
        function obj = Elliptic_GSVD(md_interface_elliptic_prior,u_in,u_out)
            obj@Randomized_GSVD(u_in,u_out);
            obj.md_interface_elliptic_prior = md_interface_elliptic_prior;
        end
   
        
    end
    
end