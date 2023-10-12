classdef Prior_Model < handle
    
    % We assume a Bayesian inverse problem with a mean zero Gaussian noise
    % model and a linear observation operator 
    
    properties
        z_prior_mean;
    end
    
     methods (Abstract, Access = public)
         
         [z_out] = Prior_Precision_Apply(this,z_in);
         
         [z_prior_mean] = Get_Prior_Mean(this);
         
     end
     
     methods (Access = public)
         
         function this = Prior_Model()
            this.z_prior_mean = this.Get_Prior_Mean();
         end
         
          function [val,grad_z] = Regularization(this,z)
            z_tmp1 = this.Prior_Precision_Apply(z-this.z_prior_mean);
            val = 0.5*(z_tmp1'*(z-this.z_prior_mean));
            grad_z = z_tmp1;
         end
         
         function [Mv] = Regularization_HessVec(this,v)
             Mv = this.Prior_Precision_Apply(v);
         end
         
         
     end
     
end

