%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef MD_Transient_Elliptic_u_Prior_Interface < MD_Scaled_u_Prior_Interface

    properties
        use_gsvd
        spatial_prior_cov
        transient_prior_cov
        u_hyperparam_interface
        determine_u_hyperparams
    end

    %% Constructor
    methods

        function this = MD_Transient_Elliptic_u_Prior_Interface(data_interface, spatial_prior_cov, transient_prior_cov)

            arguments
                data_interface MD_Data_Interface
                spatial_prior_cov MD_u_Prior_Interface
                transient_prior_cov MD_Transient_Prior_Covariance
            end
            this@MD_Scaled_u_Prior_Interface(transient_prior_cov.u_hyperparam_interface.alpha_u);

            this.use_gsvd = false;
            try
                if ~isempty(spatial_prior_cov.sing_vals)
                    this.use_gsvd = true;
                end
            end

            this.spatial_prior_cov = spatial_prior_cov;
            this.transient_prior_cov = transient_prior_cov;
            this.u_hyperparam_interface = transient_prior_cov.u_hyperparam_interface;
            this.determine_u_hyperparams = MD_Determine_u_Hyperparameters(data_interface, this.u_hyperparam_interface);

            if this.u_hyperparam_interface.adapt_time_variance
                this.determine_u_hyperparams.Determine_alpha_t(this);
            end
            this.transient_prior_cov.Set_alpha_t(this.u_hyperparam_interface.alpha_t);

            if this.u_hyperparam_interface.alpha_u == 0.0
                this.determine_u_hyperparams.Determine_alpha_u(this);
            end
            this.Set_alpha_u(this.u_hyperparam_interface.alpha_u);
        end

    end

    %% Implementation of base class virtual functions
    methods

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

            if this.use_gsvd

                u_out = 0.0 * u_in;
                for k = 1:size(u_out, 2)
                    u_tmp = reshape(u_in(:, k), this.transient_prior_cov.n_y, this.transient_prior_cov.n_t);
                    u_tmp = this.spatial_prior_cov.sing_vecs_output' * u_tmp * this.transient_prior_cov.evecs;
                    aleph = (this.spatial_prior_cov.sing_vals.^2) * this.transient_prior_cov.evals';
                    aleph = aleph ./ (1 + scalar * aleph);
                    u_tmp = u_tmp .* aleph;
                    u_tmp = this.spatial_prior_cov.sing_vecs_output * u_tmp * this.transient_prior_cov.evecs';
                    u_out(:, k) = u_tmp(:);
                end

            else

                % This version of the function avoids explicit use of the spatial prior GSVD
                u_out = 0 * u_in;
                for k = 1:size(u_in, 2)
                    u_tmp1 = reshape(u_in(:, k), this.transient_prior_cov.n_y, this.transient_prior_cov.n_t);
                    u_tmp2 = u_tmp1 * this.transient_prior_cov.evecs * diag(sqrt(this.transient_prior_cov.evals));
                    u_tmp3 = 0 * u_tmp2;
                    for j = 1:this.transient_prior_cov.n_t
                        u_tmp3(:, j) = this.spatial_prior_cov.Apply_W_u_Acute_Plus_scalar_M_u_Inverse(u_tmp2(:, j), scalar * this.transient_prior_cov.evals(j));
                    end
                    u_tmp4 = u_tmp3 * diag(sqrt(this.transient_prior_cov.evals)) * this.transient_prior_cov.evecs';
                    u_out(:, k) = u_tmp4(:);
                end

            end
        end

        function [u_out] = Apply_W_u_Acute_Inverse(this, u_in)

            if this.use_gsvd

                u_tmp = reshape(u_in, this.transient_prior_cov.n_y, this.transient_prior_cov.n_t, size(u_in, 2));
                tmp1 = pagemtimes(this.spatial_prior_cov.sing_vecs_output * diag(this.spatial_prior_cov.sing_vals.^2) * this.spatial_prior_cov.sing_vecs_output', u_tmp);
                tmp2 = pagemtimes(this.transient_prior_cov.evecs * diag(this.transient_prior_cov.evals) * this.transient_prior_cov.evecs', 'transpose', tmp1, 'transpose');
                tmp3 = pagetranspose(tmp2);
                u_out = reshape(tmp3, size(u_in, 1), size(u_in, 2));

            else

                u_out = 0.0 * u_in;
                for k = 1:size(u_out, 2)
                    u_tmp = reshape(u_in(:, k), this.transient_prior_cov.n_y, this.transient_prior_cov.n_t);
                    u_tmp = u_tmp * this.transient_prior_cov.evecs * diag(this.transient_prior_cov.evals) * this.transient_prior_cov.evecs';
                    u_tmp = this.spatial_prior_cov.Apply_W_u_Acute_Inverse(u_tmp);
                    u_out(:, k) = u_tmp(:);
                end

            end
        end

        function [u_out] = Sample_with_Covariance_W_u_Acute_Inverse(this, num_samples)

            if this.use_gsvd

                r_s = length(this.spatial_prior_cov.sing_vals);
                r_t = length(this.transient_prior_cov.evals);
                m = size(this.spatial_prior_cov.sing_vecs_output, 1) * size(this.transient_prior_cov.evecs, 1);

                omega = randn(r_s, r_t, num_samples);
                tmp1 = pagemtimes(this.spatial_prior_cov.sing_vecs_output * diag(this.spatial_prior_cov.sing_vals), omega);
                tmp2 = pagemtimes(this.transient_prior_cov.evecs * diag(sqrt(this.transient_prior_cov.evals)), 'none', tmp1, 'transpose');
                tmp3 = pagetranspose(tmp2);
                u_out = reshape(tmp3, m, num_samples);

                % This original code loops over the columns of u_out, and
                % consequently, is an easier piece of code to read.
                % The optimized code above computes the same matrix multiply
                % using tensor operations for greater efficiency
                % u_out = 0.0 * zeros(m, num_samples);
                % for k = 1:size(u_out, 2)
                %     omega = diag(this.spatial_prior_cov.sing_vals) * randn(r_s, r_t) * diag(sqrt(this.transient_prior_cov.evals));
                %     u_space = this.spatial_prior_cov.sing_vecs_output * omega;
                %     u_tmp = u_space * this.transient_prior_cov.evecs';
                %     u_out(:, k) = u_tmp(:);
                % end

            else

                % This version of the function avoids explicit use of the spatial prior GSVD
                n_t = this.transient_prior_cov.n_t;
                spatial_samples = this.spatial_prior_cov.Sample_with_Covariance_W_u_Acute_Inverse(num_samples * n_t);
                m = n_t * this.transient_prior_cov.n_y;
                u_out = 0.0 * zeros(m, num_samples);
                for k = 1:size(u_out, 2)
                    I = ((k - 1) * n_t + 1):(k * n_t);
                    u_tmp = spatial_samples(:, I) * diag(sqrt(this.transient_prior_cov.evals)) * this.transient_prior_cov.evecs';
                    u_out(:, k) = u_tmp(:);
                end

            end
        end

        function [u_out] = Sample_with_Covariance_W_u_Acute_Plus_scalar_M_u_Inverse(this, num_samples, scalar)

            if this.use_gsvd

                aleph = (this.spatial_prior_cov.sing_vals.^2) * this.transient_prior_cov.evals';
                aleph = aleph ./ (1 + scalar * aleph);

                r_s = length(this.spatial_prior_cov.sing_vals);
                r_t = length(this.transient_prior_cov.evals);
                m = size(this.spatial_prior_cov.sing_vecs_output, 1) * size(this.transient_prior_cov.evecs, 1);
                u_out = 0.0 * zeros(m, num_samples);
                for k = 1:size(u_out, 2)
                    omega = sqrt(aleph) .* randn(r_s, r_t);
                    u_space = this.spatial_prior_cov.sing_vecs_output * omega;
                    u_tmp = u_space * this.transient_prior_cov.evecs';
                    u_out(:, k) = u_tmp(:);
                end

            else

                % This version of the function avoids explicit use of the spatial prior GSVD

                n_y = this.transient_prior_cov.n_y;
                n_t = this.transient_prior_cov.n_t;
                u_out = zeros(n_t * n_y, num_samples);

                spatial_samples = zeros(n_y, n_t, num_samples);
                for j = 1:n_t
                    spatial_samples(:, j, :) = this.spatial_prior_cov.Sample_with_Covariance_W_u_Acute_Plus_scalar_M_u_Inverse(num_samples, scalar * this.transient_prior_cov.evals(j));
                end

                for k = 1:num_samples
                    u_tmp = spatial_samples(:, :, k) * diag(sqrt(this.transient_prior_cov.evals)) * this.transient_prior_cov.evecs';
                    u_out(:, k) = u_tmp(:);
                end

            end
        end

    end

end
