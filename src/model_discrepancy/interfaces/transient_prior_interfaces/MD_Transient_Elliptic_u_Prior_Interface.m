classdef MD_Transient_Elliptic_u_Prior_Interface < MD_Scaled_u_Prior_Interface

    properties
        spatial_prior_cov
        transient_prior_cov
    end

    methods

        % spatial_prior_cov should be an object of type MD_Elliptic_u_Prior_Interface
        % transient_prior_cov should be an object of type MD_Transient_Prior_Covariance
        function this = MD_Transient_Elliptic_u_Prior_Interface(spatial_prior_cov, transient_prior_cov)
            this@MD_Scaled_u_Prior_Interface(transient_prior_cov.hyperparams.alpha_u);
            this.spatial_prior_cov = spatial_prior_cov;
            this.transient_prior_cov = transient_prior_cov;
            if transient_prior_cov.hyperparams.alpha_u == 0.0
                transient_prior_cov.hyperparams.Determine_alpha_u(this);
            end
            this.Set_alpha_u(transient_prior_cov.hyperparams.alpha_u);
        end

        function [u_out] = Apply_M_u(this, u_in)
            if size(u_in, 1) == this.transient_prior_cov.n_y
                u_out = this.spatial_prior_cov.Apply_M_u(u_in);
            else
                u_out = 0.0 * u_in;
                for k = 1:size(u_out, 2)
                    u_tmp = reshape(u_in(:, k), this.transient_prior_cov.n_y, this.transient_prior_cov.n_t);
                    u_tmp = u_tmp * this.transient_prior_cov.M_t;
                    u_tmp = this.spatial_prior_cov.Apply_M_u(u_tmp);
                    u_out(:, k) = u_tmp(:);
                end
            end
        end

        function [u_out] = Apply_W_u_Acute_Plus_scalar_M_u_Inverse(this, u_in, scalar)
            u_out = 0.0 * u_in;
            for k = 1:size(u_out, 2)
                u_tmp = reshape(u_in(:, k), this.transient_prior_cov.n_y, this.transient_prior_cov.n_t);
                u_tmp = this.spatial_prior_cov.sing_vecs_output' * u_tmp * this.transient_prior_cov.M_t_inv_evecs;
                aleph = (this.spatial_prior_cov.sing_vals.^2) * this.transient_prior_cov.evals';
                aleph = aleph ./ (1 + scalar * aleph);
                u_tmp = u_tmp .* aleph;
                u_tmp = this.spatial_prior_cov.sing_vecs_output * u_tmp * this.transient_prior_cov.M_t_inv_evecs';
                u_out(:, k) = u_tmp(:);
            end
        end

        function [u_out] = Apply_W_u_Acute_Inverse(this, u_in)
            u_out = 0.0 * u_in;
            for k = 1:size(u_out, 2)
                u_tmp = reshape(u_in(:, k), this.transient_prior_cov.n_y, this.transient_prior_cov.n_t);
                u_tmp = u_tmp * this.transient_prior_cov.M_t_inv_evecs * diag(this.transient_prior_cov.evals) * this.transient_prior_cov.M_t_inv_evecs';
                u_tmp = this.spatial_prior_cov.sing_vecs_output * diag(this.spatial_prior_cov.sing_vals.^2) * this.spatial_prior_cov.sing_vecs_output' * u_tmp;
                u_out(:, k) = u_tmp(:);
            end
        end

        % Compute samples from a mean zero Gaussian with covariance W_u^{-1}
        function [u_out] = Sample_with_Covariance_W_u_Inverse(this, num_samples)
            r_s = length(this.spatial_prior_cov.sing_vals);
            r_t = length(this.transient_prior_cov.evals);
            m = size(this.spatial_prior_cov.sing_vecs_output, 1) * size(this.transient_prior_cov.evecs, 1);
            u_out = 0.0 * zeros(m, num_samples);
            for k = 1:size(u_out, 2)
                omega = diag(this.spatial_prior_cov.sing_vals) * randn(r_s, r_t) * diag(sqrt(this.transient_prior_cov.evals));
                u_space = sqrt(this.alpha_u) * this.spatial_prior_cov.sing_vecs_output * omega;
                u_tmp = u_space * this.transient_prior_cov.M_t_inv_evecs';
                u_out(:, k) = u_tmp(:);
            end
        end

        % Compute samples from a mean zero Gaussian with covariance (W_u+scalar*M_u)^{-1}
        function [u_out] = Sample_with_Covariance_W_u_Plus_scalar_M_u_Inverse(this, num_samples, scalar)
            aleph = (this.spatial_prior_cov.sing_vals.^2) * this.transient_prior_cov.evals';
            aleph = aleph ./ (1 + this.alpha_u * scalar * aleph);

            r_s = length(this.spatial_prior_cov.sing_vals);
            r_t = length(this.transient_prior_cov.evals);
            m = size(this.spatial_prior_cov.sing_vecs_output, 1) * size(this.transient_prior_cov.evecs, 1);
            u_out = 0.0 * zeros(m, num_samples);
            for k = 1:size(u_out, 2)
                omega = sqrt(aleph) .* randn(r_s, r_t);
                u_space = sqrt(this.alpha_u) * this.spatial_prior_cov.sing_vecs_output * omega;
                u_tmp = u_space * this.transient_prior_cov.M_t_inv_evecs';
                u_out(:, k) = u_tmp(:);
            end
        end

    end

end
