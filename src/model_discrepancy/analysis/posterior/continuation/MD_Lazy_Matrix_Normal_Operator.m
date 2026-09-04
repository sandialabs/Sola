%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef MD_Lazy_Matrix_Normal_Operator < handle

    % Lazy sampler for
    %
    %   Y ~ MN_{m,k}(0, W_u^{-1}, I_k).
    %
    % This object supports consistent repeated evaluations of
    %
    %   v -> Y v,
    %   u -> Y^T u,
    %
    % without explicitly materializing Y.
    %
    % Internally, it maintains
    %
    %   Q_l' W_u Q_l = I,
    %   Q_r' Q_r     = I,
    %
    % and a coupling matrix K.  Once both bases are complete, the sampled
    % matrix is
    %
    %   Y = Q_l K Q_r'.
    %
    % While the bases are incomplete, this representation is the restriction
    % of one consistent matrix-normal sample to all queried directions.

    properties
        u_prior_interface

        output_dim
        input_dim

        Q_l      % output_dim x num_left, W_u-orthonormal
        Q_r      % input_dim x num_right, Euclidean-orthonormal
        K        % num_left x num_right

        tol
    end

    methods

        function this = MD_Lazy_Matrix_Normal_Operator(u_prior_interface, input_dim, output_dim, tol)
            arguments
                u_prior_interface MD_u_Prior_Interface
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
            this.input_dim = input_dim;
            this.output_dim = output_dim;
            this.tol = tol;

            this.Q_l = zeros(output_dim, 0);
            this.Q_r = zeros(input_dim, 0);
            this.K = zeros(0, 0);
        end

        function y = Forward_Apply(this, v)
            % Compute y = Y v.

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
            % Compute y = Y^T u.

            u = u(:);

            if length(u) ~= this.output_dim
                error('Adjoint_Apply input has wrong dimension.');
            end

            % The left basis is W_u-orthonormal.  Enriching with
            %
            %   query = W_u^{-1} u
            %
            % makes the coefficient vector equal to
            %
            %   c = Q_l' u.
            query = this.u_prior_interface.Apply_W_u_Inverse(u);

            [c, added] = this.Forward_Enrich_Left(query);

            if added
                this.Backward_Enrich_Right();
            end

            y = this.Q_r * (this.K' * c);
        end

        function [c, added] = Forward_Enrich_Right(this, x)
            % Enrich the right/input basis using Euclidean orthonormality.
            %
            % Because the right covariance is I_k, the right precision is also
            % I_k and no covariance transform is needed for a forward query.

            x = x(:);

            if length(x) ~= this.input_dim
                error('Forward_Enrich_Right input has wrong dimension.');
            end

            if isempty(this.Q_r)
                c = zeros(0, 1);
                r = x;
            else
                c = this.Q_r' * x;
                r = x - this.Q_r * c;

                % One reorthogonalization pass for numerical robustness.
                dc = this.Q_r' * r;
                c = c + dc;
                r = r - this.Q_r * dc;
            end

            norm_r = norm(r);
            added = norm_r > this.tol;

            if added
                q_new = r / norm_r;
                this.Q_r = [this.Q_r, q_new];
                c = [c; norm_r];

                % New right basis vector means a new column of K.  Its
                % couplings against already-existing left directions are iid
                % standard normal variables.
                n_left = size(this.K, 1);
                new_col = randn(n_left, 1);
                this.K = [this.K, new_col];
            end
        end

        function [c, added] = Forward_Enrich_Left(this, x)
            % Enrich the left/output basis using W_u-orthonormality.
            %
            % x should be in covariance-transformed coordinates.  In the adjoint
            % apply, x = W_u^{-1} u, so
            %
            %   c = Q_l' W_u x = Q_l' u.

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

                % One W_u-weighted reorthogonalization pass for robustness.
                Wr = this.u_prior_interface.Apply_W_u(r);
                dc = this.Q_l' * Wr;
                c = c + dc;
                r = r - this.Q_l * dc;
            end

            norm_r = this.W_Norm(r);
            added = norm_r > this.tol;

            if added
                q_new = r / norm_r;
                this.Q_l = [this.Q_l, q_new];
                c = [c; norm_r];

                % New left basis vector means a new row of K.  Its couplings
                % against already-existing right directions are iid standard
                % normal variables.
                n_right = size(this.K, 2);
                new_row = randn(1, n_right);
                this.K = [this.K; new_row];
            end
        end

        function added = Backward_Enrich_Left(this)
            % Enrich the left/output side after discovering a new right
            % direction.
            %
            % Draw s ~ N(0, W_u^{-1}), remove its currently represented
            % W_u-orthogonal projection, and append the residual direction.
            %
            % The appended row of K is zero except for the bottom-right entry,
            % which stores the W_u norm of the sampled residual.

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

                % One W_u-weighted reorthogonalization pass for robustness.
                Wr = this.u_prior_interface.Apply_W_u(r);
                dc = this.Q_l' * Wr;
                r = r - this.Q_l * dc;
            end

            norm_r = this.W_Norm(r);
            added = norm_r > this.tol;

            if added
                q_new = r / norm_r;
                this.Q_l = [this.Q_l, q_new];

                n_right = size(this.K, 2);
                new_row = zeros(1, n_right);
                this.K = [this.K; new_row];

                if n_right > 0
                    this.K(end, end) = norm_r;
                end
            end
        end

        function added = Backward_Enrich_Right(this)
            % Enrich the right/input side after discovering a new left
            % direction.
            %
            % Since the right covariance is identity, draw a standard normal
            % vector and Euclidean-orthogonalize it against Q_r.
            %
            % The appended column of K is zero except for the bottom-right
            % entry, which stores the Euclidean norm of the sampled residual.

            s = randn(this.input_dim, 1);

            if isempty(this.Q_r)
                r = s;
            else
                r = s - this.Q_r * (this.Q_r' * s);

                % One Euclidean reorthogonalization pass for robustness.
                dc = this.Q_r' * r;
                r = r - this.Q_r * dc;
            end

            norm_r = norm(r);
            added = norm_r > this.tol;

            if added
                q_new = r / norm_r;
                this.Q_r = [this.Q_r, q_new];

                n_left = size(this.K, 1);
                new_col = zeros(n_left, 1);
                this.K = [this.K, new_col];

                if n_left > 0
                    this.K(end, end) = norm_r;
                end
            end
        end

        function X_cache = Explicit_Matrix(this)
            % Return the currently represented dense matrix
            %
            %   X_cache = Q_l K Q_r'.
            %
            % If both bases are complete, this is the full sampled matrix.
            % Otherwise, this is the restriction of the sampled operator to the
            % explored subspaces.

            X_cache = this.Q_l * this.K * this.Q_r';
        end

        function [kl, kr] = Basis_Dimensions(this)
            kl = size(this.Q_l, 2);
            kr = size(this.Q_r, 2);
        end

        function val = W_Inner(this, x, y)
            % Weighted inner product x' W_u y.

            val = x' * this.u_prior_interface.Apply_W_u(y);
        end

        function val = W_Norm(this, x)
            % Weighted norm sqrt(x' W_u x), protected against small negative
            % roundoff.

            inner_val = this.W_Inner(x, x);
            inner_val = real(inner_val);
            val = sqrt(max(inner_val, 0));
        end

    end

end