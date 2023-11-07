classdef Transient_Inv_Prob_Likelihood_Model < Likelihood_Model

    properties
        sigma
        obs_vec
        space_time_obs_vec
        m
        N
    end

    methods (Access = public)

        function [d_out] = Noise_Precision_Apply(this, d_in)
            d_out = diag(1 / this.sigma^2) * d_in;
        end

        function [d_out] = Observation_Operator_Apply(this, u_in)
            d_out = u_in(this.space_time_obs_vec, :);
        end

        function [u_out] = Observation_Operator_Transpose_Apply(this, d_in)
            u_out = zeros(this.m * this.N, size(d_in, 2));
            u_out(this.space_time_obs_vec, :) = d_in;
        end

        function [d] = Get_Observed_Data(this)
            u_data = load('Obs_Data.mat', 'u_data').u_data;
            d = u_data(this.space_time_obs_vec);
        end

    end

    methods (Access = public)

        function this = Transient_Inv_Prob_Likelihood_Model(m, N)
            this.sigma = 5.0;
            this.obs_vec = round(linspace(1, m, 15));
            this.space_time_obs_vec = zeros(length(this.obs_vec), N);
            for k = 1:N
                this.space_time_obs_vec(:, k) = this.obs_vec + (k - 1) * m;
            end
            this.space_time_obs_vec = this.space_time_obs_vec(:);
            this.m = m;
            this.N = N;
        end

    end

end
