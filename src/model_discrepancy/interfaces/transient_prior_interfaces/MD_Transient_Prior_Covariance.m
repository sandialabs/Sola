classdef MD_Transient_Prior_Covariance < handle

    properties
        u_hyperparam_interface
        determine_u_hyperparams
        alpha_t
        beta_t
        M_t
        S_t
        E_t
        evecs
        evals
        n_t
        n_y
    end

    methods

        function [] = Set_alpha_t(this, alpha_t_new)
            this.alpha_t = alpha_t_new;
            this.Compute_Time_Covariance_GEVP();
        end

        function [] = Set_beta_t(this, beta_t_new)
            this.beta_t = beta_t_new;
            this.E_t = this.beta_t * this.S_t + this.M_t;
        end

        function this = MD_Transient_Prior_Covariance(M_t, S_t, n_y, data_interface, u_hyperparam_interface)
            this.u_hyperparam_interface = u_hyperparam_interface;
            this.determine_u_hyperparams = MD_Determine_u_Hyperparameters(data_interface,u_hyperparam_interface);
            this.M_t = M_t;
            this.S_t = S_t;
            this.n_t = size(M_t, 1);
            this.n_y = n_y;

            if this.u_hyperparam_interface.beta_t == 0.0
                this.determine_u_hyperparams.Determine_beta_t();
            end
            this.Set_beta_t(this.u_hyperparam_interface.beta_t);
        end

        function [] = Compute_Time_Covariance_GEVP(this)
            E_t_inv = diag(sqrt(this.alpha_t)) * linsolve(this.E_t, diag(sqrt(this.alpha_t)) * eye(this.n_t));
            A = this.M_t * E_t_inv * this.M_t;
            [V, Lambda] = eig(A, this.M_t);
            n = sqrt(diag(V' * this.M_t * V));

            this.evecs = V * diag(1 ./ n);
            this.evals = diag(Lambda);
        end

        function [samples] = Sample_Time_Series(this, num_samples)
            samples = this.evecs * diag(sqrt(this.evals)) * randn(this.n_t, num_samples);
        end

    end

end
