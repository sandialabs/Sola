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
        relative_tol
        max_basis_vector_norm
        local_orthogonality_tol
        max_reorthogonalization_passes

        % Sample completion (see "Posterior sample completion" in the paper).
        % After completion, the caches are frozen and
        %   X = S + Y Q_r' Sigma + Q_l T' (I - Q_r Q_r' Sigma),
        % with S = S_left * S_right' stored through thin factors.
        auto_complete = 0   % 0: off; f>0: complete once (k_f+k_t)(n_u+r) >= f*n_u*r (f=1 is break-even)
        is_completed = false
        S_left     % output_dim x num_comp, P_u S_u, P_u = I - Q_l Q_l' W_u
        S_right    % input_dim x num_comp, G with G G' = Sigma - Sigma Q_r Q_r' Sigma
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
            this.relative_tol = sqrt(tol);
            this.max_basis_vector_norm = 1.e6;
            this.local_orthogonality_tol = 1.e-6;
            this.max_reorthogonalization_passes = 3;

            this.Q_l = zeros(output_dim, 0);
            this.W_Q_l = zeros(output_dim, 0);
            this.Q_r = zeros(input_dim, 0);
            this.Sigma_Q_r = zeros(input_dim, 0);
            this.Y = zeros(output_dim, 0);
            this.T = zeros(input_dim, 0);
        end


        function y = Forward_Apply(this, v)

            v = v(:);
            Sigma_v = this.sigma_apply(v);

            if this.is_completed
                y = this.Completed_Forward_Apply(v, Sigma_v);
                return;
            end

            total_norm = sqrt(max(real(v' * Sigma_v), 0));

            if isempty(this.Q_r)
                c = zeros(0, 1);
                r = v;
                Sigma_r = Sigma_v;
            else
                c = this.Q_r' * Sigma_v;
                r = v - this.Q_r * c;
                Sigma_r = Sigma_v - this.Sigma_Q_r * c;
                [r, Sigma_r, c] = this.Orthogonalize_Right_Residual(r, Sigma_r, c, total_norm);
            end

            norm_r = sqrt(max(real(r' * Sigma_r), 0));
            add_basis = this.Should_Add_Basis_Vector(r, norm_r, total_norm);

            if add_basis
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

            this.Check_Auto_Complete();

        end


        function y = Adjoint_Apply(this, u)
            % Compute y = X' u.
            u = u(:);

            if this.is_completed
                y = this.Completed_Adjoint_Apply(u);
                return;
            end

            x = this.u_prior_interface.Apply_W_u_Inverse(u);
            total_norm = sqrt(max(real(x' * u), 0));

            if isempty(this.Q_l)
                c = zeros(0, 1);
                r = x;
                W_r = u;
            else
                c = this.Q_l' * u;
                r = x - this.Q_l * c;
                W_r = u - this.W_Q_l * c;
                [r, W_r, c] = this.Orthogonalize_Left_Residual(r, W_r, c, total_norm);
            end

            norm_r = sqrt(max(real(r' * W_r), 0));
            add_basis = this.Should_Add_Basis_Vector(r, norm_r, total_norm);

            if add_basis
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

            this.Check_Auto_Complete();

        end

        % ------------------------------------------------------------
        % Sample completion
        % ------------------------------------------------------------

        function Check_Auto_Complete(this)
            if this.auto_complete > 0 && ~this.is_completed
                [kl, kr] = this.Basis_Dimensions();
                if (kl + kr) * (this.input_dim + this.output_dim) >= this.auto_complete * this.input_dim * this.output_dim
                    this.Complete_Sample();
                end
            end
        end

        function Complete_Sample(this)
            % Draw the part of X not determined by the caches, conditionally
            % on everything revealed so far, and freeze the caches.
            if this.is_completed
                return;
            end

            % Factor G with G G' = Sigma - Sigma Q_r Q_r' Sigma. Sigma is
            % formed through input_dim applies and eigendecomposed.
            I_r = eye(this.input_dim);
            Sigma = zeros(this.input_dim);
            for j = 1:this.input_dim
                Sigma(:, j) = this.sigma_apply(I_r(:, j));
            end
            Sigma = (Sigma + Sigma') / 2;
            Sigma_res = Sigma - this.Sigma_Q_r * this.Sigma_Q_r';
            Sigma_res = (Sigma_res + Sigma_res') / 2;
            [U, L] = eig(Sigma_res, 'vector');
            keep = L > max([L; 0]) * eps * this.input_dim * 1e2;
            G = U(:, keep) * diag(sqrt(L(keep)));
            num_comp = size(G, 2);

            % P_u S_u with S_u columns ~ N(0, W_u^{-1})
            if num_comp > 0
                S_u = this.u_prior_interface.Sample_with_Covariance_W_u_Inverse(num_comp);
                if ~isempty(this.Q_l)
                    S_u = S_u - this.Q_l * (this.W_Q_l' * S_u);
                    S_u = S_u - this.Q_l * (this.W_Q_l' * S_u); % re-projection
                end
            else
                S_u = zeros(this.output_dim, 0);
            end
            this.S_left = S_u;
            this.S_right = G;
            this.is_completed = true;
        end

        function y = Completed_Forward_Apply(this, v, Sigma_v)
            % y = S v + Y Q_r' Sigma v + Q_l T' (v - Q_r Q_r' Sigma v)
            y = this.S_left * (this.S_right' * v);
            if ~isempty(this.Q_r)
                c = this.Q_r' * Sigma_v;
                y = y + this.Y * c;
                if ~isempty(this.Q_l)
                    y = y + this.Q_l * (this.T' * (v - this.Q_r * c));
                end
            elseif ~isempty(this.Q_l)
                y = y + this.Q_l * (this.T' * v);
            end
        end

        function y = Completed_Adjoint_Apply(this, u)
            % y = S' u + Sigma Q_r Y' u + (I - Sigma Q_r Q_r') T Q_l' u
            y = this.S_right * (this.S_left' * u);
            if ~isempty(this.Q_r)
                y = y + this.Sigma_Q_r * (this.Y' * u);
            end
            if ~isempty(this.Q_l)
                w = this.T * (this.Q_l' * u);
                if ~isempty(this.Q_r)
                    w = w - this.Sigma_Q_r * (this.Q_r' * w);
                end
                y = y + w;
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
            for jj = 1:this.max_reorthogonalization_passes
                deta = this.W_Q_l' * eta;
                if norm(deta) <= this.tol
                    break;
                end
                eta = eta - this.Q_l * deta;
            end

        end

        function eta = Sample_Right_Residual(this)
            s = this.sigma_sample(1);
            if isempty(this.Q_r)
                eta = s;
                return;
            end
            eta = s - this.Sigma_Q_r * (this.Q_r' * s);
            for jj = 1:this.max_reorthogonalization_passes
                deta = this.Q_r' * eta;
                if norm(deta) <= this.tol
                    break;
                end
                eta = eta - this.Sigma_Q_r * deta;
            end
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

        function metrics = Cache_Health_Metrics(this)
            metrics.left_orth_err = 0.0;
            metrics.right_orth_err = 0.0;
            metrics.cross_err = 0.0;

            if ~isempty(this.Q_l)
                metrics.left_orth_err = norm(this.Q_l' * this.W_Q_l - eye(size(this.Q_l, 2)), 'fro');
            end

            if ~isempty(this.Q_r)
                metrics.right_orth_err = norm(this.Q_r' * this.Sigma_Q_r - eye(size(this.Q_r, 2)), 'fro');
            end

            if ~isempty(this.Q_l) && ~isempty(this.Q_r)
                left_cache = this.W_Q_l' * this.Y;
                right_cache = this.T' * this.Q_r;
                metrics.cross_err = norm(left_cache - right_cache, 'fro') / ...
                    max([1, norm(left_cache, 'fro'), norm(right_cache, 'fro')]);
            end
        end

        function tf = Should_Add_Basis_Vector(this, r, norm_r, total_norm)
            tf = norm_r > max([this.tol, this.relative_tol * total_norm, norm(r) / this.max_basis_vector_norm]);
        end

        function [r, Sigma_r, c] = Orthogonalize_Right_Residual(this, r, Sigma_r, c, total_norm)
            threshold = max(this.tol, this.relative_tol * total_norm);
            for jj = 1:this.max_reorthogonalization_passes
                norm_r = sqrt(max(real(r' * Sigma_r), 0));
                dc = this.Q_r' * Sigma_r;
                if norm_r <= threshold || norm(dc) <= norm_r * this.local_orthogonality_tol
                    return;
                end
                c = c + dc;
                r = r - this.Q_r * dc;
                Sigma_r = Sigma_r - this.Sigma_Q_r * dc;
            end
        end

        function [r, W_r, c] = Orthogonalize_Left_Residual(this, r, W_r, c, total_norm)
            threshold = max(this.tol, this.relative_tol * total_norm);
            for jj = 1:this.max_reorthogonalization_passes
                norm_r = sqrt(max(real(r' * W_r), 0));
                dc = this.W_Q_l' * r;
                if norm_r <= threshold || norm(dc) <= norm_r * this.local_orthogonality_tol
                    return;
                end
                c = c + dc;
                r = r - this.Q_l * dc;
                W_r = W_r - this.W_Q_l * dc;
            end
        end

    end

end
