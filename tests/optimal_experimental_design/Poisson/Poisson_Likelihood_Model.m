classdef Poisson_Likelihood_Model < Likelihood_Model

    properties
        sigma
        obs_vec
        m
    end

    methods (Access = public)

        function [d_out] = Noise_Precision_Apply(this, d_in)
            d_out = diag(1 / this.sigma^2) * d_in;
        end

        function [d_out] = Observation_Operator_Apply(this, u_in)
            d_out = u_in(this.obs_vec, :);
        end

        function [u_out] = Observation_Operator_Transpose_Apply(this, d_in)
            u_out = zeros(this.m, size(d_in, 2));
            u_out(this.obs_vec, :) = d_in;
        end

        function [d] = Get_Observed_Data(this)
            d = zeros(length(this.obs_vec), 1);
        end

    end

    methods (Access = public)

        function this = Poisson_Likelihood_Model(sigma, obs_vec, m)
            this.sigma = sigma;
            this.obs_vec = obs_vec;
            this.m = m;
        end

    end

end
