classdef Likelihood_Model < handle

    % We assume a Bayesian inverse problem with a mean zero Gaussian noise
    % model and a linear observation operator

    properties

    end

    methods (Abstract, Access = public)

        [d_out] = Noise_Precision_Apply(this, d_in)

        [d_out] = Observation_Operator_Apply(this, u_in)

        [u_out] = Observation_Operator_Transpose_Apply(this, d_in)

        [d] = Get_Observed_Data(this)

    end

    methods (Access = public)

        function this = Likelihood_Model()

        end

        function [val, grad_u] = Misfit(this, u)
            u_tmp1 = this.Observation_Operator_Apply(u) - this.Get_Observed_Data();
            u_tmp2 = this.Noise_Precision_Apply(u_tmp1);
            val = 0.5 * (u_tmp1' * u_tmp2);
            grad_u = this.Observation_Operator_Transpose_Apply(u_tmp2);
        end

        function [Mv] = Misfit_HessVec(this, v)
            u_tmp1 = this.Observation_Operator_Apply(v);
            u_tmp2 = this.Noise_Precision_Apply(u_tmp1);
            Mv = this.Observation_Operator_Transpose_Apply(u_tmp2);
        end

    end

end
