classdef HDSA_Sabl_Transient_Prior_Covariance < HDSA_Transient_Prior_Covariance

    properties

    end

    methods

        function this = HDSA_Sabl_Transient_Prior_Covariance(beta_tu, beta_iu, beta_td, T, N)
            h = T / (N - 1);

            M_t = diag(4 * ones(1, N)) + diag(ones(1, N - 1), 1) + diag(ones(1, N - 1), -1);
            M_t(1, 1) = .5 * M_t(1, 1);
            M_t(end, end) = .5 * M_t(end, end);
            M_t = (1 / 6) * h * M_t;

            S_t = diag(2 * ones(1, N)) + (-1) * diag(ones(1, N - 1), 1) + (-1) * diag(ones(1, N - 1), -1);
            S_t(1, 1) = .5 * S_t(1, 1);
            S_t(end, end) = .5 * S_t(end, end);
            S_t = (1 / h) * S_t;

            this@HDSA_Transient_Prior_Covariance(beta_tu, beta_iu, beta_td, M_t, S_t);
            num_evals = N;
            oversampling = 20;
            this.Compute_Time_Covariance_GEVP(num_evals, oversampling);
        end

    end

end
