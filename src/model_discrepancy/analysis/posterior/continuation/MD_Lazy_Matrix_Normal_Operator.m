%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef MD_Lazy_Matrix_Normal_Operator < handle

    % Symmetric lazy sampler for
    %
    %   X ~ MN_{m,r}(0, W_u^{-1}, Sigma_r),
    %
    % where Sigma_r is available through matrix-vector products and sampling.
    %
    % This implementation uses the symmetric revealed-action cache
    %
    %   Y = X * Q_r,
    %   T = X' * W_u * Q_l,
    %

    properties
        u_prior_interface
        sigma_apply
        sigma_sample

        output_dim
        input_dim

        Q_l        % output_dim x num_left, W_u-orthonormal
        W_Q_l
        Q_r        % input_dim x num_right, Sigma_r-orthonormal
        Sigma_Q_r  % input_dim x num_right, Sigma_r * Q_r
        Y          % output_dim x num_right, X * Q_r
        T          % input_dim x num_left, X' * W_u * Q_l

        tol
    end

    methods

        function this = MD_Lazy_Matrix_Normal_Operator(u_prior_interface, sigma_apply, sigma_sample, input_dim, output_dim, tol)
            arguments
                u_prior_interface MD_u_Prior_Interface
                sigma_apply
                sigma_sample
                input_dim (1, 1) {mustBeNumeric}
                output_dim (1, 1) {mustBeNumeric}
                tol (1, 1) {mustBeNumeric} = 1e-10
            end

            if input_dim < 0 || floor(input_dim) ~= input_dim
                error('input_dim must be a nonnegative integer.');
            end

            if output_dim < 0 || floor(output_dim) ~= output_dim
                error('output_dim must be a nonnegative integer.');
            end

            this.u_prior_interface = u_prior_interface;
            this.sigma_apply = sigma_apply;
            this.sigma_sample = sigma_sample;
            this.input_dim = input_dim;
            this.output_dim = output_dim;
            this.tol = tol;

            this.Q_l = zeros(output_dim, 0);
            this.W_Q_l = zeros(output_dim, 0);
            this.Q_r = zeros(input_dim, 0);
            this.Sigma_Q_r = zeros(input_dim, 0);
            this.Y = zeros(output_dim, 0);
            this.T = zeros(input_dim, 0);
        end


        function y = Forward_Apply(this, v)

            Sigma_v = this.sigma_apply(v);

            if isempty(this.Q_r)
                c = zeros(0, 1);
                r = v;
                Sigma_r = Sigma_v;
            else
                c = this.Q_r' * Sigma_v;
                r = v - this.Q_r * c;
                Sigma_r = Sigma_v - this.Sigma_Q_r * c;

                % One Sigma_r-weighted reorthogonalization pass.
                dc = this.Q_r' * Sigma_r;
                c = c + dc;
                r = r - this.Q_r * dc;
                Sigma_r = Sigma_r - this.Sigma_Q_r * dc;
            end

            norm_r = sqrt(max(real(r' * Sigma_r), 0));

            if norm_r > this.tol
                q_new = r / norm_r;
                Sigma_q_new = Sigma_r / norm_r;

                if isempty(this.Q_l)
                    y_cond = zeros(this.output_dim, 1);
                else
                    left_coeffs = this.T' * q_new;
                    y_cond = this.Q_l * left_coeffs;
                end

                y_new = y_cond + this.Sample_Left_Residual();
                this.Q_r = [this.Q_r, q_new];
                this.Sigma_Q_r = [this.Sigma_Q_r, Sigma_q_new];
                this.Y = [this.Y, y_new];
                c = [c; norm_r];

            end

            if isempty(this.Y)
                y = zeros(this.output_dim, 1);
            else
                y = this.Y * c;
            end

        end


        function y = Adjoint_Apply(this, u)
            % Compute y = X' u.
            x = this.u_prior_interface.Apply_W_u_Inverse(u);

            if isempty(this.Q_l)
                c = zeros(0, 1);
                r = x;
                W_r = u;
            else
                c = this.Q_l' * u;
                r = x - this.Q_l * c;
                W_r = u - this.W_Q_l * c;

                % One W_u-weighted reorthogonalization pass.
                dc = this.W_Q_l' * r;
                c = c + dc;
                r = r - this.Q_l * dc;
                W_r = W_r - this.W_Q_l * dc;
            end

            norm_r = sqrt(max(real(r' * W_r), 0));

            if norm_r > this.tol
                q_new = r / norm_r;
                W_q_new = W_r / norm_r;

                if isempty(this.Q_r)
                    t_cond = zeros(this.input_dim, 1);
                else
                    right_coeffs = this.Y' * W_q_new;
                    t_cond = this.Sigma_Q_r * right_coeffs;
                end

                t_new = t_cond + this.Sample_Right_Residual();
                this.Q_l = [this.Q_l, q_new];
                this.W_Q_l = [this.W_Q_l, W_q_new];
                this.T = [this.T, t_new];
                c = [c; norm_r];
            end

            if isempty(this.T)
                y = zeros(this.input_dim, 1);
            else
                y = this.T * c;
            end

        end

        % ------------------------------------------------------------
        % Residual samplers
        % ------------------------------------------------------------

        function eta = Sample_Left_Residual(this)
            s = this.u_prior_interface.Sample_with_Covariance_W_u_Inverse(1);
            if isempty(this.Q_l)
                eta = s;
                return;
            end
            eta = s - this.Q_l * (this.W_Q_l' * s);
            % One reorthogonalization pass.
            eta = eta - this.Q_l * (this.W_Q_l' * eta);

        end

        function eta = Sample_Right_Residual(this)
            s = this.sigma_sample(1);
            if isempty(this.Q_r)
                eta = s;
                return;
            end
            eta = s - this.Sigma_Q_r * (this.Q_r' * s);
            % One reorthogonalization pass.
            eta = eta - this.Sigma_Q_r * (this.Q_r' * eta);
        end

        function [kl, kr] = Basis_Dimensions(this)
            kl = size(this.Q_l, 2);
            kr = size(this.Q_r, 2);
        end

        function val = W_Inner(this, x, y)
            val = x' * this.u_prior_interface.Apply_W_u(y);
        end

        function val = W_Norm(this, x)
            val = sqrt(max(real(this.W_Inner(x, x)), 0));
        end
    end

end
