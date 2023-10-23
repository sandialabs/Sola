classdef Thermal_Prior_Model < Inf_Dim_Prior_Model
    
    
    properties
        con;
        L;
    end
    
    methods (Access = public)
        
        function [z_out] = Laplacian_Like_Apply(this,z_in)
            z_out = this.L*z_in;
        end
        
        function [z_out] = Laplacian_Like_Transpose_Apply(this,z_in)
            z_out = this.L'*z_in;
        end
        
        function [z_out] = Laplacian_Like_Inverse_Apply(this,z_in)
            z_out = linsolve(this.L,z_in);
        end
        
        function [z_out] = Laplacian_Like_Transpose_Inverse_Apply(this,z_in)
            z_out = linsolve(this.L',z_in);
        end

        function [z_out] = Mass_Matrix_Apply(this,z_in)
            z_out = this.con.M*z_in;
        end
        
        function [z_out] = Mass_Matrix_Inverse_Apply(this,z_in)
            z_out = linsolve(this.con.M,z_in);
        end
        
        function [z_prior_mean] = Get_Prior_Mean(this)
            z_prior_mean = 2*ones(this.con.m,1);
        end
    end
    
    methods (Access = public)
        
        function this = Thermal_Prior_Model(con)
           this.con = con; 
           this.L = (.5)*( (4.e-2)*con.S + con.M );
        end
        
        
    end
    
end

