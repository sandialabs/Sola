classdef Adv_Diff_Likelihood_Model < Likelihood_Model

    properties
        noise_precision
        obs_operator
    end

    methods (Access = public)

        function [d_out] = Noise_Precision_Apply(this, d_in)
            d_out = this.noise_precision * d_in;
        end

        function [d_out] = Observation_Operator_Apply(this, u_in)
            d_out = this.obs_operator * u_in;
        end

        function [u_out] = Observation_Operator_Transpose_Apply(this, d_in)
            u_out = this.obs_operator' * d_in;
        end

        function [d] = Get_Observed_Data(this)
            u_data = load('Obs_Data.mat', 'u_data').u_data;
            d = this.obs_operator * u_data;
        end

    end

    methods (Access = public)

        function this = Adv_Diff_Likelihood_Model(m)
            obs_locations = round(linspace(1, m, 8)');
            obs_operator = eye(m);
            obs_operator = obs_operator(obs_locations, :);
            this.noise_precision = (1.e6) * (1.0 / 3.0)^2 * eye(length(obs_locations));
            this.obs_operator = obs_operator;
        end

    end

end
