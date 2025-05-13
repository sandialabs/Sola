classdef MD_Elliptic_u_Prior_Interface < MD_Scaled_u_Prior_Interface

    %%%%%%%%%%%%%%%%% Defined covariance as a squared inverse elliptic operator %%%%%%%%%%%%%%%%%

    properties
        sing_vecs_output
        sing_vals
    end

    %% Pure virtual functions for user implementation
    methods (Abstract, Access = public)

        [u_out] = Apply_E_u_Inverse(this, u_in)

        [u_out] = Apply_E_u_Inverse_Transpose(this, u_in)

        [u_out] = Apply_M_u(this, u_in)

    end

    %% Constructor and helper functions
    methods

        function this = MD_Elliptic_u_Prior_Interface(alpha_u)
            arguments
                alpha_u (1, 1) {mustBeNumeric}
            end
            this@MD_Scaled_u_Prior_Interface(alpha_u);
        end

        function [] = Compute_E_u_Inverse_GSVD(this, num_sing_vals, oversampling, num_subspace_iters, u_vec)
            gsvd = Elliptic_GSVD(this, u_vec, u_vec);
            [~, this.sing_vecs_output, this.sing_vals] = gsvd.Compute_GSVD(num_sing_vals, oversampling, num_subspace_iters);
        end

    end

    %% Implementation of base class functions
    methods

        function [u_out] = Apply_W_u_Acute_Plus_scalar_M_u_Inverse(this, u_in, scalar)
            K = (this.sing_vals.^2) ./ (1 + scalar * this.sing_vals.^2);
            u_out = this.sing_vecs_output * diag(K) * this.sing_vecs_output' * u_in;
        end

        function [u_out] = Apply_W_u_Acute_Inverse(this, u_in)
            u_out = this.sing_vecs_output * diag(this.sing_vals.^2) * this.sing_vecs_output' * u_in;
        end

        function [u_out] = Sample_with_Covariance_W_u_Acute_Inverse(this, num_samples)
            r = length(this.sing_vals);
            u_out = this.sing_vecs_output * diag(this.sing_vals) * randn(r, num_samples);
        end

        function [u_out] = Sample_with_Covariance_W_u_Acute_Plus_scalar_M_u_Inverse(this, num_samples, scalar)
            K = (this.sing_vals.^2) ./ (1 + scalar * this.sing_vals.^2);
            r = length(this.sing_vals);
            u_out = this.sing_vecs_output * diag(sqrt(K)) * randn(r, num_samples);
        end

    end

end
