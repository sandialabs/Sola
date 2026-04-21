%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef MD_Transient_Prior_Covariance < handle

    properties
        u_hyperparam_interface
        determine_u_hyperparams
        alpha_t
        beta_t
        M_t
        S_t
        E_t
        W_t
        evecs
        evals
        n_t
        n_y
    end

    %% Constructor and helper functions
    methods

        function this = MD_Transient_Prior_Covariance(M_t, S_t, n_y, data_interface, u_hyperparam_interface)

            arguments
                M_t (:, :) {mustBeNumeric}
                S_t (:, :) {mustBeNumeric}
                n_y (1, 1) {mustBeNumeric}
                data_interface MD_Data_Interface
                u_hyperparam_interface MD_u_Hyperparameter_Interface
            end
            this.u_hyperparam_interface = u_hyperparam_interface;
            this.determine_u_hyperparams = MD_Determine_u_Hyperparameters(data_interface, u_hyperparam_interface);
            this.M_t = M_t;
            this.S_t = S_t;
            this.n_t = size(M_t, 1);
            this.n_y = n_y;

            if this.u_hyperparam_interface.beta_t == 0.0
                this.determine_u_hyperparams.Determine_beta_t();
            end
            this.Set_beta_t(this.u_hyperparam_interface.beta_t);
        end

        function [] = Set_alpha_t(this, alpha_t_new)
            this.alpha_t = alpha_t_new;
            this.Compute_Time_Covariance_GEVP();
        end

        function [] = Set_beta_t(this, beta_t_new)
            this.beta_t = beta_t_new;
            this.E_t = this.beta_t * this.S_t + this.M_t;
        end

        function [] = Compute_Time_Covariance_GEVP(this)
            this.W_t = diag(sqrt(1 ./ this.alpha_t)) * this.E_t * diag(sqrt(1 ./ this.alpha_t));
            [V, Lambda] = eig(this.W_t, this.M_t, 'vector');
            [~, I] = sort(Lambda, 'ascend');
            Lambda = Lambda(I);
            V = V(:, I);
            n = sqrt(diag(V' * this.M_t * V));
            this.evecs = V * diag(1 ./ n) * diag(sign(V(1, :)));
            this.evals = 1 ./ Lambda;
        end

    end

end
