classdef HDSA_Transient_Prior_Covariance < handle

    properties
        beta_tu
        beta_iu
        beta_td
        E_tu
        E_td
        evecs
        evals
    end

    methods

        function this = HDSA_Transient_Prior_Covariance(beta_tu, beta_iu, beta_td, M_t, S_t)
            this.beta_tu = beta_tu;
            this.beta_iu = beta_iu;
            this.beta_td = beta_td;
            this.E_tu = beta_tu * S_t + M_t;
            this.E_tu(1, 1) = this.E_tu(1, 1) + beta_iu;
            this.E_td = beta_td * S_t + M_t;
        end

        function [] = Compute_Time_Covariance_GEVP(this, num_evals, oversampling)
            gevp = Time_Covariance_GEVP(this.E_tu, this.E_td);
            [this.evecs, this.evals] = gevp.Compute_GEVP(num_evals, oversampling);
        end

    end

end
