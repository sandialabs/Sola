%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Poisson_Likelihood_Model < Likelihood_Model

    properties
        d
    end

    methods (Access = public)

        function [d_out] = Noise_Precision_Apply(this, d_in)
            d_out = d_in;
        end

        function [d_out] = Observation_Operator_Apply(this, u_in)
            d_out = u_in;
        end

        function [u_out] = Observation_Operator_Transpose_Apply(this, d_in)
            u_out = d_in;
        end

        function [d] = Get_Observed_Data(this)
            d = this.d;
        end

    end

    methods (Access = public)

        function this = Poisson_Likelihood_Model(con)
            this.d = con.State_Solve(2 * con.diff_coeff * ones(con.m, 1));
        end

    end

end
