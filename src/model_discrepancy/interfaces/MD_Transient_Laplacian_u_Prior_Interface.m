classdef MD_Transient_Laplacian_u_Prior_Interface < MD_Elliptic_u_Prior_Interface

    properties
        beta_u
        S_s
        M_s
        E_s
        hyperparams
        transient_prior_cov
    end

    methods (Access = public)

        function [u_out] = Apply_E_u_Inverse(this, u_in)
            u_out = linsolve(this.E_s, u_in);
        end

        function [u_out] = Apply_E_u_Inverse_Transpose(this, u_in)
            u_out = linsolve(this.E_s', u_in);
        end

        function [u_out] = Apply_M_u(this, u_in)
            if size(u_in, 1) == this.transient_prior_cov.n_y
                u_out = this.M_s * u_in;
            else
                u_out = 0.0 * u_in;
                for k = 1:size(u_out, 2)
                    u_tmp = reshape(u_in(:, k), this.transient_prior_cov.n_y, this.transient_prior_cov.n_t);
                    u_tmp = u_tmp * this.transient_prior_cov.M_t;
                    u_tmp = this.M_s * u_tmp;
                    u_out(:, k) = u_tmp(:);
                end
            end
        end

    end

    methods

        function [] = Set_beta_u(this,beta_u_new)
            this.beta_u = beta_u_new;
            this.E_s = this.beta_u * this.S_s + this.M_s;
        end

        function [] = Compute_E_s_Inverse_GSVD(this, num_sing_vals, oversampling, num_subspace_iters, u_vec)
            gsvd = Elliptic_GSVD(this, u_vec, u_vec);
            [~, this.sing_vecs_output, this.sing_vals] = gsvd.Compute_GSVD(num_sing_vals, oversampling, num_subspace_iters);
        end

        function this = MD_Transient_Laplacian_u_Prior_Interface(S,M,hyperparams,transient_prior_cov)
            this@MD_Elliptic_u_Prior_Interface(hyperparams.alpha_u)
            this.S_s = S;
            this.M_s = M;
            this.hyperparams = hyperparams;
            this.transient_prior_cov = transient_prior_cov;

            if this.hyperparams.beta_u == 0.0
                this.hyperparams.Determine_beta_u();
            end
            this.Set_beta_u(this.hyperparams.beta_u);

            m = size(this.S_s,1);
            u_vec = zeros(m,1);
            if this.hyperparams.gsvd_num_sing_vals == 0
                this.hyperparams.Determine_GSVD_Hyperparameters();
            end
            this.Compute_E_s_Inverse_GSVD(this.hyperparams.gsvd_num_sing_vals, this.hyperparams.gsvd_oversampling, this.hyperparams.gsvd_num_subspace_iter, u_vec);

            if this.hyperparams.alpha_u == 0.0
                this.hyperparams.Determine_alpha_u(this);
            end
            this.Set_alpha_u(this.hyperparams.alpha_u);

        end
    end

end
