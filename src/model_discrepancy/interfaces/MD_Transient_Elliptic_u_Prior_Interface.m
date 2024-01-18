classdef MD_Transient_Elliptic_u_Prior_Interface < MD_u_Prior_Interface

    properties
        sing_vecs_input
        sing_vecs_output
        sing_vals
        alpha_u
        transient_prior_cov
    end

    methods (Abstract, Access = public)

        %% Pure virtual functions

        [u_out] = Apply_E_u_Inverse(this, u_in)

        [u_out] = Apply_E_u_Inverse_Transpose(this, u_in)

        [u_out] = Apply_Spatial_M_u(this, u_in)

    end

    methods

        function [] = Compute_E_u_Inverse_GSVD(this, num_sing_vals, oversampling, num_subspace_iters, u_vec)
            gsvd = Elliptic_GSVD(this, u_vec, u_vec);
            [this.sing_vecs_input, this.sing_vecs_output, this.sing_vals] = gsvd.Compute_GSVD(num_sing_vals, oversampling, num_subspace_iters);
        end

        function this = MD_Transient_Elliptic_u_Prior_Interface(alpha_u, transient_prior_cov)
            this.alpha_u = alpha_u;
            this.transient_prior_cov = transient_prior_cov;
        end

        function [u_out] = Apply_M_u(this, u_in)
            if size(u_in, 1) == this.transient_prior_cov.n_y
                u_out = this.Apply_Spatial_M_u(u_in);
            else
                u_out = 0.0 * u_in;
                for k = 1:size(u_out, 2)
                    u_tmp = reshape(u_in(:, k), this.transient_prior_cov.n_y, this.transient_prior_cov.n_t);
                    u_tmp = u_tmp * this.transient_prior_cov.M_t;
                    u_tmp = this.Apply_Spatial_M_u(u_tmp);
                    u_out(:, k) = u_tmp(:);
                end
            end
        end

        function [u_out] = Apply_W_u_Plus_scalar_M_u_Inverse(this, u_in, scalar)
            u_out = 0.0 * u_in;
            for k = 1:size(u_out, 2)
                u_tmp = reshape(u_in(:, k), this.transient_prior_cov.n_y, this.transient_prior_cov.n_t);
                u_tmp = this.sing_vecs_output' * u_tmp * this.transient_prior_cov.M_t_inv_evecs;
                aleph = (this.sing_vals.^2) * this.transient_prior_cov.evals';
                aleph = aleph ./ (1 + this.alpha_u * scalar * aleph);
                u_tmp = u_tmp .* aleph;
                u_tmp = this.sing_vecs_output * u_tmp * this.transient_prior_cov.M_t_inv_evecs';
                u_out(:, k) = this.alpha_u * u_tmp(:);
            end
        end

        function [u_out] = Apply_W_u_Inverse(this, u_in)
            u_out = 0.0 * u_in;
            for k = 1:size(u_out, 2)
                u_tmp = reshape(u_in(:, k), this.transient_prior_cov.n_y, this.transient_prior_cov.n_t);
                u_tmp = u_tmp * this.transient_prior_cov.M_t_inv_evecs * diag(this.transient_prior_cov.evals) * this.transient_prior_cov.M_t_inv_evecs';
                u_tmp = this.alpha_u * this.sing_vecs_output * diag(this.sing_vals.^2) * this.sing_vecs_output' * u_tmp;
                u_out(:, k) = u_tmp(:);
            end
        end

        % Factorize W_u^{-1}=F*F^T, function gives u_out=F*u_in
        % This function must be implemented to enable posterior update sampling
        function [u_out] = Apply_W_u_Inverse_Factor(this, u_in)
            r = length(this.sing_vals);
            u_out = 0.0 * u_in;
            for k = 1:size(u_out, 2)
                u_tmp = reshape(u_in(:, k), this.transient_prior_cov.n_y, this.transient_prior_cov.n_t);
                u_tmp = u_tmp * diag(sqrt(this.transient_prior_cov.evals)) * this.transient_prior_cov.M_t_inv_evecs';
                u_tmp = sqrt(this.alpha_u) * this.sing_vecs_output * diag(this.sing_vals) * u_tmp(1:r, :);
                u_out(:, k) = u_tmp(:);
            end
        end

        % Factorize (W_u+scalar*M_u)^{-1}=F*F^T, function gives u_out=F*u_in
        % This function must be implemented to enable posterior update sampling
        function [u_out] = Apply_W_u_Plus_scalar_M_u_Inverse_Factor(this, u_in, scalar)
            aleph = (this.sing_vals.^2) * this.transient_prior_cov.evals';
            aleph = aleph ./ (1 + this.alpha_u * scalar * aleph);
            u_out = 0.0 * u_in;
            r = length(this.sing_vals);
            for k = 1:size(u_out, 2)
                u_tmp = reshape(u_in(:, k), this.transient_prior_cov.n_y, this.transient_prior_cov.n_t);
                u_tmp = u_tmp(1:r, :) .* sqrt(aleph);
                u_tmp = this.sing_vecs_output * u_tmp * this.transient_prior_cov.M_t_inv_evecs';
                u_out(:, k) = sqrt(this.alpha_u) * u_tmp(:);
            end
        end

    end

end
