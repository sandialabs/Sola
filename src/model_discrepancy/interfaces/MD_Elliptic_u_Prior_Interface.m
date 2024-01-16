classdef MD_Elliptic_u_Prior_Interface < MD_u_Prior_Interface

    properties
        sing_vecs_input
        sing_vecs_output
        sing_vals
        alpha_u
    end

    methods (Abstract, Access = public)

        %% Pure virtual functions

        [u_out] = Apply_E_u_Inverse(this, u_in)

        [u_out] = Apply_E_u_Inverse_Transpose(this, u_in)

        [u_out] = Apply_M_u(this, u_in)

    end

    methods

        function [] = Compute_E_u_Inverse_GSVD(this, num_sing_vals, oversampling, num_subspace_iters, u_vec)
            gsvd = Elliptic_GSVD(this, u_vec, u_vec);
            [this.sing_vecs_input, this.sing_vecs_output, this.sing_vals] = gsvd.Compute_GSVD(num_sing_vals, oversampling, num_subspace_iters);
        end

        function this = MD_Elliptic_u_Prior_Interface(alpha_u)
            this.alpha_u = alpha_u;
        end

        function [u_out] = Apply_W_u_Plus_scalar_M_u_Inverse(this, u_in, scalar)
            K = (this.sing_vals.^2) ./ (1 + this.alpha_u * scalar * this.sing_vals.^2);
            u_out = this.alpha_u * this.sing_vecs_output * diag(K) * this.sing_vecs_output' * u_in;
        end

        function [u_out] = Apply_W_u_Inverse(this, u_in)
            u_out = this.alpha_u * this.sing_vecs_output * diag(this.sing_vals.^2) * this.sing_vecs_output' * u_in;
        end

        % Factorize W_u^{-1}=F*F^T, function gives u_out=F*u_in
        % This function must be implemented to enable posterior update sampling
        function [u_out] = Apply_W_u_Inverse_Factor(this, u_in)
            r = length(this.sing_vals);
            u_out = sqrt(this.alpha_u) * this.sing_vecs_output * diag(this.sing_vals) * u_in(1:r, :);
        end

        % Factorize (W_u+scalar*M_u)^{-1}=F*F^T, function gives u_out=F*u_in
        % This function must be implemented to enable posterior update sampling
        function [u_out] = Apply_W_u_Plus_scalar_M_u_Inverse_Factor(this, u_in, scalar)
            K = (this.sing_vals.^2) ./ (1 + this.alpha_u * scalar * this.sing_vals.^2);
            r = length(this.sing_vals);
            u_out = sqrt(this.alpha_u) * this.sing_vecs_output * diag(sqrt(K)) * u_in(1:r, :);
        end

    end

end
