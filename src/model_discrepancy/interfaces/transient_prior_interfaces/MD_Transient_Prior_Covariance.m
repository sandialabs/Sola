classdef MD_Transient_Prior_Covariance < handle

    properties
        hyperparams
        beta_t
        M_t
        S_t
        E_t
        evecs
        evals
        M_t_inv_evecs
        n_t
        n_y
    end

    methods

        function [] = Set_beta_t(this, beta_t_new)
            this.beta_t = beta_t_new;
            this.E_t = this.beta_t * this.S_t + this.M_t;
        end

        function this = MD_Transient_Prior_Covariance(M_t, S_t, n_y, hyperparams)
            this.hyperparams = hyperparams;
            this.M_t = M_t;
            this.S_t = S_t;
            this.n_t = size(M_t, 1);
            this.n_y = n_y;

            if this.hyperparams.beta_t == 0.0
                this.hyperparams.Determine_beta_t();
            end
            this.Set_beta_t(this.hyperparams.beta_t);

            num_evals = this.n_t;
            oversampling = 0;
            this.Compute_Time_Covariance_GEVP(num_evals, oversampling);
        end

        function [] = Compute_Time_Covariance_GEVP(this, num_evals, oversampling)
            gevp = Time_Covariance_GEVP(this.E_t, this.M_t);
            [this.evecs, this.evals] = gevp.Compute_GEVP(num_evals, oversampling);
            this.M_t_inv_evecs = this.M_t \ this.evecs;
        end

        function [samples] = Sample_Time_Series(this, num_samples)
            samples = this.M_t_inv_evecs * diag(sqrt(this.evals)) * randn(this.n_t, num_samples);
        end

    end

end
