%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Forward_Operator_GSVD < Randomized_GSVD

    properties
        prior
        likelihood
        con
    end

    methods (Access = public)

        function [vec_out] = Apply_Operator(this, vec_in)
            tmp1 = -this.con.c_z_Apply(vec_in);
            tmp2 = this.con.c_u_Inverse_Apply(tmp1);
            vec_out = this.likelihood.Observation_Operator_Apply(tmp2);
        end

        function [vec_out] = Apply_Operator_Transpose(this, vec_in)
            tmp1 = this.likelihood.Observation_Operator_Transpose_Apply(vec_in);
            tmp2 = -this.con.c_u_Transpose_Inverse_Apply(tmp1);
            vec_out = this.con.c_z_Transpose_Apply(tmp2);
        end

        function [vec_out] = Apply_Input_Weighting_Operator_Inverse(this, vec_in)
            tmp1 = this.prior.Laplacian_Like_Transpose_Inverse_Apply(vec_in);
            tmp2 = this.prior.Mass_Matrix_Apply(tmp1);
            vec_out = this.prior.Laplacian_Like_Inverse_Apply(tmp2);
        end

        function [vec_out] = Apply_Output_Weighting_Operator(this, vec_in)
            vec_out = vec_in;
        end

    end

    methods

        function this = Forward_Operator_GSVD(vec_in, vec_out, prior, likelihood, con)
            this@Randomized_GSVD(vec_in, vec_out);
            this.prior = prior;
            this.likelihood = likelihood;
            this.con = con;
        end

    end

end
