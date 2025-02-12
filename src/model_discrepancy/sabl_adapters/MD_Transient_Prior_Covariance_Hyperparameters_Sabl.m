classdef MD_Transient_Prior_Covariance_Hyperparameters_Sabl < MD_Transient_Prior_Covariance_Hyperparameters

    properties

    end

    methods

        function this = MD_Transient_Prior_Covariance_Hyperparameters_Sabl(T, n_t, n_y,hyperparams)
            h = T / (n_t - 1);

            M_t = diag(4 * ones(1, n_t)) + diag(ones(1, n_t - 1), 1) + diag(ones(1, n_t - 1), -1);
            M_t(1, 1) = .5 * M_t(1, 1);
            M_t(end, end) = .5 * M_t(end, end);
            M_t = (1 / 6) * h * M_t;

            S_t = diag(2 * ones(1, n_t)) + (-1) * diag(ones(1, n_t - 1), 1) + (-1) * diag(ones(1, n_t - 1), -1);
            S_t(1, 1) = .5 * S_t(1, 1);
            S_t(end, end) = .5 * S_t(end, end);
            S_t = (1 / h) * S_t;

            this@MD_Transient_Prior_Covariance_Hyperparameters(M_t, S_t, n_y, hyperparams);
            num_evals = n_t;
            oversampling = 0;
            this.Compute_Time_Covariance_GEVP(num_evals, oversampling);
        end

    end

end