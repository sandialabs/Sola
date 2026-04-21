%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Likelihood_Model_Example_1 < Likelihood_Model

    properties

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
            d = load('Observed_Data.mat', 'd').d;
        end

        function this = Likelihood_Model_Example_1()

        end

    end

end
