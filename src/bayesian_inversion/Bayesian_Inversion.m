classdef Bayesian_Inversion < handle
    
    % We assume a Bayesian inverse problem with a mean zero Gaussian noise
    % model, a Gaussian prior, and a linear observation operator 
    
    properties
        likelihood;
        prior
        obj;
        con;
        opt;
    end
     
     methods (Access = public)
         
         function this = Bayesian_Inversion(likelihood,prior,con)
             this.likelihood = likelihood;
             this.prior = prior;
             this.obj = Bayesian_Inversion_Objective(likelihood,prior);
             this.con = con;
             this.opt = Reduced_Space_Optimization(this.obj,this.con);
         end
         
         function [u_map,z_map] = Compute_MAP_Point(this,z0)
             [u_map,z_map] = this.opt.Optimize(z0);
         end
         
     end
     
end

