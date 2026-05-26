%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Poisson_Prior_Model < Prior_Model

    properties
        poisson_con
        L
    end

    methods (Access = public)

        function this = Poisson_Prior_Model(poisson_con, norm_scale, grad_scale)
            this.poisson_con = poisson_con;
            this.L = grad_scale * poisson_con.S + norm_scale *  poisson_con.M;
        end

        function [z_out] = Prior_Covariance_Factor_Apply(this, z)
            z_out = linsolve(this.L, this.poisson_con.M * z);
        end

        function [z_out] = Prior_Precision_Apply(this, z_in)
            z_out = linsolve(this.poisson_con.M, this.L * z_in);
            z_out = linsolve(this.poisson_con.M, this.L * z_out);
        end

        function [z_out] = Prior_Covariance_Apply(this, z_in)
            temp = this.Prior_Covariance_Factor_Apply(z_in);
            z_out = this.Prior_Covariance_Factor_Apply(temp);
        end

        function [z_prior_mean] = Get_Prior_Mean(this)
            z_prior_mean = zeros(this.poisson_con.m, 1);
        end

        function [z_out] = Mass_Matrix_Apply(this, z_in)
            z_out = this.poisson_con.M * z_in;
        end

        function [z_out] = Mass_Matrix_Inverse_Apply(this, z_in)
            z_out = linsolve(this.poisson_con.M, z_in);
        end

    end

end
