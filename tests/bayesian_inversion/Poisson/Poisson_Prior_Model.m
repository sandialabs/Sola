%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Poisson_Prior_Model < Inf_Dim_Prior_Model

    properties
        poisson_con
        L
    end

    methods (Access = public)

        function [z_out] = Laplacian_Like_Apply(this, z_in)
            z_out = this.L * z_in;
        end

        function [z_out] = Laplacian_Like_Transpose_Apply(this, z_in)
            z_out = this.L' * z_in;
        end

        function [z_out] = Laplacian_Like_Inverse_Apply(this, z_in)
            z_out = linsolve(this.L, z_in);
        end

        function [z_out] = Laplacian_Like_Transpose_Inverse_Apply(this, z_in)
            z_out = linsolve(this.L', z_in);
        end

        function [z_out] = Mass_Matrix_Apply(this, z_in)
            z_out = this.poisson_con.M * z_in;
        end

        function [z_out] = Mass_Matrix_Inverse_Apply(this, z_in)
            z_out = linsolve(this.poisson_con.M, z_in);
        end

        function [z_prior_mean] = Get_Prior_Mean(this)
            z_prior_mean = 2 * this.poisson_con.diff_coeff * ones(this.poisson_con.m, 1);
        end

    end

    methods (Access = public)

        function this = Poisson_Prior_Model(poisson_con)
            this.poisson_con = poisson_con;
            this.L = (1.e-3) * poisson_con.S + poisson_con.M;
        end

    end

end
