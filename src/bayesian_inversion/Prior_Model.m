classdef Prior_Model < handle

    % We assume a Bayesian inverse problem with a mean zero Gaussian noise
    % model and a linear observation operator

    properties

    end

    methods (Abstract, Access = public)

        [z_out] = Prior_Precision_Apply(this, z_in)

        [z_out] = Prior_Covariance_Apply(this, z_in)

        [z_prior_mean] = Get_Prior_Mean(this)

        % Assume that the prior covariance is factorized as
        % Gamma = L*L^T, this function computes z_out = L*z_in
        [z_out] = Prior_Covariance_Factor_Apply(this, z_in)

    end

    methods (Access = public)

        function this = Prior_Model()

        end

        function [val, grad_z] = Regularization(this, z)
            tmp1 = z - this.Get_Prior_Mean();
            tmp2 = this.Prior_Precision_Apply(tmp1);
            val = 0.5 * (tmp2' * tmp1);
            grad_z = tmp2;
        end

        function [Mv] = Regularization_HessVec(this, v)
            Mv = this.Prior_Precision_Apply(v);
        end

        function [Z_prior] = Compute_Prior_Samples(this, num_samps)
            z_prior_mean = this.Get_Prior_Mean();
            Omega = randn(length(z_prior_mean), num_samps);
            Z_prior = this.Prior_Covariance_Factor_Apply(Omega) + z_prior_mean;
        end

    end

end
