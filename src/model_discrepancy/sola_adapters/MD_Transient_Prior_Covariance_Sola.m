%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef MD_Transient_Prior_Covariance_Sola < MD_Transient_Prior_Covariance

    properties

    end

    methods

        function this = MD_Transient_Prior_Covariance_Sola(data_interface, u_hyperparam_interface, T, n_t, n_y)
            arguments
                data_interface MD_Data_Interface
                u_hyperparam_interface MD_u_Hyperparameter_Interface
                T (1, 1) double
                n_t (1, 1) double
                n_y (1, 1) double
            end
            h = T / (n_t - 1);

            M_t = diag(4 * ones(1, n_t)) + diag(ones(1, n_t - 1), 1) + diag(ones(1, n_t - 1), -1);
            M_t(1, 1) = .5 * M_t(1, 1);
            M_t(end, end) = .5 * M_t(end, end);
            M_t = (1 / 6) * h * M_t;

            S_t = diag(2 * ones(1, n_t)) + (-1) * diag(ones(1, n_t - 1), 1) + (-1) * diag(ones(1, n_t - 1), -1);
            S_t(1, 1) = .5 * S_t(1, 1);
            S_t(end, end) = .5 * S_t(end, end);
            S_t = (1 / h) * S_t;

            this@MD_Transient_Prior_Covariance(M_t, S_t, n_y, data_interface, u_hyperparam_interface);
        end

    end

end
