classdef Hessian_GEVP < Randomized_GEVP
    
    properties
        z_opt;
        md_interface;
        evals;
        evecs;
        is_computed;
        normalization_coeff;
    end
    
    methods (Access = public)
        
        function [vec_out] = Apply_Operator(obj,vec_in)
            vec_out = obj.md_interface.Apply_RS_Hessian(vec_in,obj.z_opt);
        end
        
        function [vec_out] = Apply_Weighting_Operator(obj,vec_in)
            vec_out = (1/obj.normalization_coeff)*obj.md_interface.Apply_W_z(vec_in);
        end
        
        function [vec_out] = Apply_Weighting_Operator_Inverse(obj,vec_in)
           vec_out = obj.normalization_coeff*obj.md_interface.Apply_W_z_Inverse(vec_in);
        end
               
        function [vec_out] = Apply_Weighting_Operator_Inverse_Factor(obj,vec_in)
            vec_out = sqrt(obj.normalization_coeff)*obj.md_interface.Apply_W_z_Inverse_Factor(vec_in);
        end
        
    end
    
    methods
        function obj = Hessian_GEVP(md_interface,z_opt)
            obj@Randomized_GEVP(z_opt);
            obj.z_opt = z_opt;
            obj.md_interface = md_interface;
            obj.is_computed = false;
            obj.normalization_coeff = md_interface.Apply_W_z(z_opt)'*z_opt;
        end
   
        function [] = Compute_Hessian_GEVP(obj,num_evals,oversampling)
            [obj.evecs,obj.evals] = obj.Compute_GEVP(num_evals, oversampling);
            obj.is_computed = true;
        end
        
       function [z_out] = Apply_Projected_RS_Hessian_Inverse(obj,z_in)
            z_out = obj.evecs*diag(1./obj.evals)*obj.evecs'*z_in;
       end
        
    end
    
end