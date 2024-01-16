classdef MD_Transient_Prior_Covariance < handle

    properties
        beta_t
        beta_i
        E_t
        M_t
        evecs
        evals
        M_t_inv_evecs
        n_t
        n_y
    end

    methods

        function this = MD_Transient_Prior_Covariance(beta_t, beta_i, M_t, S_t, n_y)
            this.beta_t = beta_t;
            this.beta_i = beta_i;
            this.E_t = beta_t * S_t + M_t;
            this.E_t(1, 1) = this.E_t(1, 1) + beta_i;
            this.M_t = M_t;
            this.n_t = size(M_t, 1);
            this.n_y = n_y;
        end

        function [] = Compute_Time_Covariance_GEVP(this, num_evals, oversampling)
            gevp = Time_Covariance_GEVP(this.E_t, this.M_t);
            [this.evecs, this.evals] = gevp.Compute_GEVP(num_evals, oversampling);
            this.M_t_inv_evecs = this.M_t \ this.evecs;
        end

    end

end
