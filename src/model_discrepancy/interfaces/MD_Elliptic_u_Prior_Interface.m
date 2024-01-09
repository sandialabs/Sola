classdef MD_Elliptic_u_Prior_Interface < MD_u_Prior_Interface

    properties
        sing_vecs_input
        sing_vecs_output
        sing_vals
        alpha_u
        is_transient
        transient_prior_cov
    end

    methods (Abstract, Access = public)

        %% Pure virtual functions

        [u_out] = Apply_E_u_Inverse(this, u_in)

        [u_out] = Apply_E_u_Inverse_Transpose(this, u_in)

        [u_out] = Apply_M_u(this, u_in)

        [u_out] = Apply_M_u_Inverse(this, u_in)

        [u_out] = Apply_E_d(this, u_in)

        [u_out] = Apply_E_d_Transpose(this, u_in)

    end

    methods

        function [] = Compute_E_u_Inverse_GSVD(this, num_sing_vals, oversampling, num_subspace_iters, u_vec)
            gsvd = Elliptic_GSVD(this, u_vec, u_vec);
            [this.sing_vecs_input, this.sing_vecs_output, this.sing_vals] = gsvd.Compute_GSVD(num_sing_vals, oversampling, num_subspace_iters);
        end

        function this = MD_Elliptic_u_Prior_Interface(alpha_u, transient_prior_cov)
            this.alpha_u = alpha_u;
            switch nargin
                case 1
                    this.is_transient = false;
                case 2
                    if isempty(transient_prior_cov)
                        this.is_transient = false;
                    else
                        this.is_transient = true;
                        this.transient_prior_cov = transient_prior_cov;
                    end
            end
        end

        function [u_out] = Apply_W_d(this, u_in)
            if ~this.is_transient
                tmp1 = this.Apply_E_d(u_in);
                tmp2 = this.Apply_M_u_Inverse(tmp1);
                u_out = this.Apply_E_d_Transpose(tmp2);
            elseif size(u_in, 1) == this.transient_prior_cov.n_y
                tmp1 = this.Apply_E_d(u_in);
                tmp2 = this.Apply_M_u_Inverse(tmp1);
                u_out = this.Apply_E_d_Transpose(tmp2);
            else
                u_tmp = reshape(u_in, this.transient_prior_cov.n_y, this.transient_prior_cov.n_t);
                u_tmp = u_tmp * this.transient_prior_cov.E_td;
                tmp1 = this.Apply_E_d(u_tmp);
                tmp2 = this.Apply_M_u_Inverse(tmp1);
                u_out = this.Apply_E_d_Transpose(tmp2);
                u_out = u_out(:);
            end
        end

        function [u_out] = Apply_W_u_Plus_scalar_W_d_Inverse(this, u_in, scalar)
            if ~this.is_transient
                K = (this.sing_vals.^2) ./ (1 + this.alpha_u * scalar * this.sing_vals.^2);
                u_out = this.alpha_u * this.sing_vecs_output * diag(K) * this.sing_vecs_output' * u_in;
            else
                u_tmp = reshape(u_in, this.transient_prior_cov.n_y, this.transient_prior_cov.n_t);
                u_tmp = this.sing_vecs_output' * u_tmp * this.transient_prior_cov.E_td_inv_evecs;
                aleph = (this.sing_vals.^2) * this.transient_prior_cov.evals';
                aleph = aleph ./ (1 + this.alpha_u * scalar * aleph);
                u_tmp = u_tmp .* aleph;
                u_tmp = this.sing_vecs_output * u_tmp * this.transient_prior_cov.E_td_inv_evecs';
                u_out = this.alpha_u * u_tmp(:);
            end
        end

        function [u_out] = Apply_W_u_Inverse(this, u_in)
            if ~this.is_transient
                u_out = this.alpha_u * this.sing_vecs_output * diag(this.sing_vals.^2) * this.sing_vecs_output' * u_in;
            else
                u_tmp = reshape(u_in, this.transient_prior_cov.n_y, this.transient_prior_cov.n_t);
                u_tmp = u_tmp * this.transient_prior_cov.E_td_inv_evecs * diag(this.transient_prior_cov.evals) * this.transient_prior_cov.E_td_inv_evecs';
                u_out = this.alpha_u * this.sing_vecs_output * diag(this.sing_vals.^2) * this.sing_vecs_output' * u_tmp;
                u_out = u_out(:);
            end
        end

        % Factorize W_u^{-1}=F*F^T, function gives u_out=F*u_in
        % This function must be implemented to enable posterior update sampling
        function [u_out] = Apply_W_u_Inverse_Factor(this, u_in)
            if ~this.is_transient
                r = length(this.sing_vals);
                u_out = sqrt(this.alpha_u) * this.sing_vecs_output * diag(this.sing_vals) * u_in(1:r, :);
            else
                r = length(this.sing_vals);
                u_out = 0.0 * u_in;
                for k = 1:size(u_out, 2)
                    u_tmp = reshape(u_in(:, k), this.transient_prior_cov.n_y, this.transient_prior_cov.n_t);
                    u_tmp = u_tmp * diag(sqrt(this.transient_prior_cov.evals)) * this.transient_prior_cov.E_td_inv_evecs';
                    u_tmp = sqrt(this.alpha_u) * this.sing_vecs_output * diag(this.sing_vals) * u_tmp(1:r, :);
                    u_out(:, k) = u_tmp(:);
                end
            end
        end

        % Factorize (W_u+scalar*W_d)^{-1}=F*F^T, function gives u_out=F*u_in
        % This function must be implemented to enable posterior update sampling
        function [u_out] = Apply_W_u_Plus_scalar_W_d_Inverse_Factor(this, u_in, scalar)
            if ~this.is_transient
                K = (this.sing_vals.^2) ./ (1 + this.alpha_u * scalar * this.sing_vals.^2);
                r = length(this.sing_vals);
                u_out = sqrt(this.alpha_u) * this.sing_vecs_output * diag(sqrt(K)) * u_in(1:r, :);
            else
                aleph = (this.sing_vals.^2) * this.transient_prior_cov.evals';
                aleph = aleph ./ (1 + this.alpha_u * scalar * aleph);
                u_out = 0.0 * u_in;
                r = length(this.sing_vals);
                for k = 1:size(u_out, 2)
                    u_tmp = reshape(u_in(:, k), this.transient_prior_cov.n_y, this.transient_prior_cov.n_t);
                    u_tmp = u_tmp(1:r, :) .* sqrt(aleph);
                    u_tmp = this.sing_vecs_output * u_tmp * this.transient_prior_cov.E_td_inv_evecs';
                    u_out(:, k) = sqrt(this.alpha_u) * u_tmp(:);
                end
            end
        end

    end

end
