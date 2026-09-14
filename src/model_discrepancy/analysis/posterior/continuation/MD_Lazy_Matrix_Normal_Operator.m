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

                this.Q_r = [this.Q_r, q_new];
                this.Sigma_Q_r = [this.Sigma_Q_r, Sigma_q_new];
                c = [c; norm_r];

                new_col = this.T' * q_new;
                this.K = [this.K, new_col];

                % %
                s = this.u_prior_interface.Sample_with_Covariance_W_u_Inverse(1);

                if isempty(this.Q_l)
                    rh = s;
                else
                    Ws = this.u_prior_interface.Apply_W_u(s);
                    rh = s - this.Q_l * (this.Q_l' * Ws);

                    Wr = this.u_prior_interface.Apply_W_u(rh);
                    dc = this.Q_l' * Wr;
                    rh = rh - this.Q_l * dc;
                end

                norm_rh = this.W_Norm(rh);
                if norm_rh > this.tol
                    q_new = rh / norm_rh;
                    k_row = zeros(1, size(this.K, 2));
                    if size(this.K, 2) > 0
                        k_row(end) = norm_rh;
                    end
                    this.Append_Left_Basis(q_new, k_row);
                end
            end

            y = this.Q_l * (this.K * c);
        end

        function y = Adjoint_Apply(this, u)
            % Compute y = X' u.
            x = this.u_prior_interface.Apply_W_u_Inverse(u);

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

            if norm_r > this.tol
                q_new = r / norm_r;
                c = [c; norm_r];
                this.Append_Left_Basis(q_new, zeros(1, size(this.K, 2)));
            end

            y = this.T * c;
        end

        function [] = Append_Left_Basis(this, q_new, k_row)

            s = this.sigma_sample(1);

            if isempty(this.Q_r)
                t_new = s;
            else
                r = s - this.Sigma_Q_r * (this.Q_r' * s);
                % Reduce roundoff error
                r = r - this.Sigma_Q_r * (this.Q_r' * r);
                t_new = this.Sigma_Q_r * k_row(:) + r;
            end

            this.Q_l = [this.Q_l, q_new];
            this.T = [this.T, t_new];
            this.K = [this.K; k_row(:)'];
        end

        function X_cache = Explicit_Matrix(this)
            X_cache = this.Q_l * this.T';
        end

        function [kl, kr] = Basis_Dimensions(this)
            kl = size(this.Q_l, 2);
            kr = size(this.Q_r, 2);
        end

        function val = W_Inner(this, x, y)
            val = x' * this.u_prior_interface.Apply_W_u(y);
        end

        function val = W_Norm(this, x)
            val = sqrt(max(this.W_Inner(x, x), 0));
        end

    end

end

