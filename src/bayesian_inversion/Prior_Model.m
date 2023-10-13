classdef Prior_Model < handle
    
    % We assume a Bayesian inverse problem with a mean zero Gaussian noise
    % model and a linear observation operator 
    
    properties

    end
    
     methods (Abstract, Access = public)
         
         [z_out] = Prior_Precision_Apply(this,z_in);
         
         [z_prior_mean] = Get_Prior_Mean(this);
         
     end
     
     methods (Access = public)
         
         function this = Prior_Model()

         end
         
         function [val,grad_z] = Regularization(this,z)
             tmp1 = z-this.Get_Prior_Mean();
             tmp2 = this.Prior_Precision_Apply(tmp1);
             val = 0.5*(tmp2'*tmp1);
             grad_z = tmp2;
         end
         
         function [Mv] = Regularization_HessVec(this,v)
             Mv = this.Prior_Precision_Apply(v);
         end
         
         
     end
     
end

