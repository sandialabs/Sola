classdef HDSA_MD_Prior_Sampling < handle
    
    properties
        md_interface;
        u_opt;
        z_opt;
    end
        
    methods
        function obj = HDSA_MD_Prior_Sampling(md_interface)
            obj.md_interface = md_interface;
            obj.u_opt = obj.md_interface.Load_Optimal_u();
            obj.z_opt = obj.md_interface.Load_Optimal_z();
        end
           
        function [z_samples] = Prior_z_Samples(obj,num_samps)
            n = length(obj.z_opt);
            Omega = randn(n,num_samps);
            z_samples = obj.md_interface.Apply_W_z_Inverse_Factor(Omega);
        end
        
        function [delta_samples] = Prior_Discrepancy_Samples_at_z_opt(obj,num_samps)
            m = length(obj.u_opt);
            Omega = randn(m,num_samps);
            delta_samples = obj.md_interface.Apply_W_u_Inverse_Factor(Omega);
        end
        
        function [delta_samples] = Prior_Discrepancy_Samples(obj,z,num_samps)
            Z = z-obj.z_opt;
            Sigma = Z'*obj.md_interface.Apply_W_z_Inverse(Z);
            m = length(obj.u_opt);
            p = size(Z,2);
            R = chol(Sigma);
            
            delta_samples = cell(num_samps,1);
            for k = 1:num_samps
                u_vec = randn(m,p)*R + randn(m,1);
                delta_samples{k} = obj.md_interface.Apply_W_u_Inverse_Factor(u_vec);
            end
        end
            
      
    end
    
end