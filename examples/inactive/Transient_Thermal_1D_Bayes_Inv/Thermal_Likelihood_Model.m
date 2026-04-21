%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Thermal_Likelihood_Model < Likelihood_Model

    properties
        sigma
        obs_vec
        space_time_obs_vec
        n_y
        n_t
    end

    methods (Access = public)

        function [d_out] = Noise_Precision_Apply(this, d_in)
            d_out = diag(1 / this.sigma^2) * d_in;
        end

        function [d_out] = Observation_Operator_Apply(this, u_in)
            d_out = u_in(this.space_time_obs_vec, :);
        end

        function [u_out] = Observation_Operator_Transpose_Apply(this, d_in)
            u_out = zeros(this.n_y * this.n_t, size(d_in, 2));
            u_out(this.space_time_obs_vec, :) = d_in;
        end

        function [d] = Get_Observed_Data(this)
            u_data = load('Obs_Data.mat', 'u_data').u_data;
            d = u_data(this.space_time_obs_vec);
        end

    end

    methods (Access = public)

        function this = Thermal_Likelihood_Model(n_y, n_t)
            this.sigma = 5.0;
            this.obs_vec = round(linspace(1, n_y, 15));
            this.space_time_obs_vec = zeros(length(this.obs_vec), n_t);
            for k = 1:n_t
                this.space_time_obs_vec(:, k) = this.obs_vec + (k - 1) * n_y;
            end
            this.space_time_obs_vec = this.space_time_obs_vec(:);
            this.n_y = n_y;
            this.n_t = n_t;
        end

    end

end
