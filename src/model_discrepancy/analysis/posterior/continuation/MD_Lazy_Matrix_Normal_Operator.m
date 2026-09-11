%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef MD_Lazy_Matrix_Normal_Operator < handle

    % Lazy sampler for
    %
    %   X ~ MN_{m,r}(0, W_u^{-1}, Sigma_r),
    %
    % where Sigma_r is available through matrix-vector products and sampling.
    % The right basis is Sigma_r-orthonormal,
    %
    %   Q_r' Sigma_r Q_r = I,
    %
    % and the represented sample is
    %
    %   X = Q_l K Q_r' Sigma_r.
    %
    % The sampler also stores T = X' W_u Q_l.  This lets adjoint evaluations
    % return beta-space vectors without completing the entire right basis.

    properties
        u_prior_interface
        sigma_apply
        sigma_sample

        output_dim
        input_dim

        Q_l        % output_dim x num_left, W_u-orthonormal
        Q_r        % input_dim x num_right, Sigma_r-orthonormal
        Sigma_Q_r  % input_dim x num_right, Sigma_r * Q_r
        T          % input_dim x num_left, X' * W_u * Q_l
        K          % num_left x num_right

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
            this.Q_r = zeros(input_dim, 0);
            this.Sigma_Q_r = zeros(input_dim, 0);
            this.T = zeros(input_dim, 0);
            this.K = zeros(0, 0);
        end

        function y = Forward_Apply(this, v)
            % Compute y = X v.

            v = v(:);

            if length(v) ~= this.input_dim
                error('Forward_Apply input has wrong dimension.');
            end

            [c, added] = this.Forward_Enrich_Right(v);

            if added
                this.Backward_Enrich_Left();
            end

            y = this.Q_l * (this.K * c);
        end

        function y = Adjoint_Apply(this, u)
            % Compute y = X' u.

            u = u(:);

            if length(u) ~= this.output_dim
                error('Adjoint_Apply input has wrong dimension.');
            end

            query = this.u_prior_interface.Apply_W_u_Inverse(u);
            c = this.Forward_Enrich_Left(query);

            y = this.T * c;
        end

        function [c, added] = Forward_Enrich_Right(this, x)
            % Enrich the right/input basis using Sigma_r-orthonormality.

            x = x(:);

            if length(x) ~= this.input_dim
                error('Forward_Enrich_Right input has wrong dimension.');
            end

            Sigma_x = this.Apply_Sigma(x);

            if isempty(this.Q_r)
                c = zeros(0, 1);
                r = x;
                Sigma_r = Sigma_x;
            else
                c = this.Q_r' * Sigma_x;
                r = x - this.Q_r * c;
                Sigma_r = Sigma_x - this.Sigma_Q_r * c;

                % One Sigma_r-weighted reorthogonalization pass.
                dc = this.Q_r' * Sigma_r;
                c = c + dc;
                r = r - this.Q_r * dc;
                Sigma_r = Sigma_r - this.Sigma_Q_r * dc;
            end

            norm_r = sqrt(max(real(r' * Sigma_r), 0));
            added = norm_r > this.tol;

            if added
                q_new = r / norm_r;
                Sigma_q_new = Sigma_r / norm_r;

                this.Q_r = [this.Q_r, q_new];
                this.Sigma_Q_r = [this.Sigma_Q_r, Sigma_q_new];
                c = [c; norm_r];

                new_col = this.T' * q_new;
                this.K = [this.K, new_col];
            end
        end

        function [c, added] = Forward_Enrich_Left(this, x)
            % Enrich the left/output basis using W_u-orthonormality.

            x = x(:);

            if length(x) ~= this.output_dim
                error('Forward_Enrich_Left input has wrong dimension.');
            end

            if isempty(this.Q_l)
                c = zeros(0, 1);
                r = x;
            else
                Wx = this.u_prior_interface.Apply_W_u(x);
                c = this.Q_l' * Wx;
                r = x - this.Q_l * c;

                Wr = this.u_prior_interface.Apply_W_u(r);
                dc = this.Q_l' * Wr;
                c = c + dc;
                r = r - this.Q_l * dc;
            end

            norm_r = this.W_Norm(r);
            added = norm_r > this.tol;

            if added
                q_new = r / norm_r;
                c = [c; norm_r];

                n_right = size(this.K, 2);
                k_row = randn(1, n_right);
                this.Append_Left_Basis(q_new, k_row);
            end
        end

        function added = Backward_Enrich_Left(this)
            % Enrich the left/output side after discovering a new right
            % direction by drawing a W_u^{-1} sample residual.

            s = this.u_prior_interface.Sample_with_Covariance_W_u_Inverse(1);
            s = s(:, 1);

            if length(s) ~= this.output_dim
                error('Sample_with_Covariance_W_u_Inverse returned a vector with wrong dimension.');
            end

            if isempty(this.Q_l)
                r = s;
            else
                Ws = this.u_prior_interface.Apply_W_u(s);
                r = s - this.Q_l * (this.Q_l' * Ws);

                Wr = this.u_prior_interface.Apply_W_u(r);
                dc = this.Q_l' * Wr;
                r = r - this.Q_l * dc;
            end

            norm_r = this.W_Norm(r);
            added = norm_r > this.tol;

            if added
                q_new = r / norm_r;

                n_right = size(this.K, 2);
                k_row = zeros(1, n_right);

                if n_right > 0
                    k_row(end) = norm_r;
                end

                this.Append_Left_Basis(q_new, k_row);
            end
        end

        function [] = Append_Left_Basis(this, q_new, k_row)

            n_right = size(this.Q_r, 2);

            if length(k_row) ~= n_right
                error('Append_Left_Basis prescribed row has wrong dimension.');
            end

            t_new = this.Sample_Right_Conditional(k_row(:));

            this.Q_l = [this.Q_l, q_new];
            this.T = [this.T, t_new];
            this.K = [this.K; k_row(:)'];
        end

        function t = Sample_Right_Conditional(this, k_col)
            % Draw t = X' W_u q_l conditioned on q_j' t = K(l,j).

            s = this.Sample_Sigma(1);
            s = s(:, 1);

            if isempty(this.Q_r)
                t = s;
            else
                r = s - this.Sigma_Q_r * (this.Q_r' * s);

                % One pass to reduce roundoff in the interpolation conditions
                % Q_r' r = 0.
                r = r - this.Sigma_Q_r * (this.Q_r' * r);
                t = this.Sigma_Q_r * k_col + r;
            end
        end

        function X_cache = Explicit_Matrix(this)
            X_cache = this.Q_l * this.T';
        end

        function [kl, kr] = Basis_Dimensions(this)
            kl = size(this.Q_l, 2);
            kr = size(this.Q_r, 2);
        end

        function Sigma_x = Apply_Sigma(this, x)
            Sigma_x = this.sigma_apply(x);
            Sigma_x = Sigma_x(:);

            if length(Sigma_x) ~= this.input_dim
                error('sigma_apply returned a vector with wrong dimension.');
            end
        end

        function samples = Sample_Sigma(this, num_samples)
            samples = this.sigma_sample(num_samples);

            if size(samples, 1) ~= this.input_dim || size(samples, 2) ~= num_samples
                error('sigma_sample returned an array with wrong dimensions.');
            end
        end

        function val = W_Inner(this, x, y)
            val = x' * this.u_prior_interface.Apply_W_u(y);
        end

        function val = W_Norm(this, x)
            val = sqrt(max(this.W_Inner(x, x), 0));
        end

    end

end
