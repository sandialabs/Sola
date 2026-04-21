%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Matrix_Sqrt < handle

    properties

    end

    methods (Abstract, Access = public)

        %% Pure virtual functions
        [vec_out] = Matrix_Apply(this, vec_in)

    end

    methods (Access = public)

        %% Preconditioner G such that G^T*G \approx A^{-1}
        % Overload these functions if a preconditioner is available

        % Apply G * vec_in
        function [vec_out] = Preconditioner_Apply(this, vec_in)
            vec_out = vec_in;
        end

        % Apply G^T * vec_in
        function [vec_out] = Preconditioner_Transpose_Apply(this, vec_in)
            vec_out = vec_in;
        end

        % Apply G^{-1} * vec_in
        function [vec_out] = Preconditioner_Inverse_Apply(this, vec_in)
            vec_out = vec_in;
        end

    end

    methods

        function this = Matrix_Sqrt()

        end

        function [vec_out, relres] = Matrix_Sqrt_Apply(this, vec_in)
            n = size(vec_in, 1);
            d = size(vec_in, 2);
            vec_out = zeros(n, d);
            A = @(v) this.Preconditioner_Apply(this.Matrix_Apply(this.Preconditioner_Transpose_Apply(v)));
            tol = 1.e-8;
            relres = cell(d, 1);
            for k = 1:d
                [tmp, relres{k}] = this.Krylov_Sqrt(A, vec_in(:, k), n, tol);
                vec_out(:, k) = this.Preconditioner_Inverse_Apply(tmp);
            end

            % Note that the preconditioner implies that
            % Matrix_Sqrt_Apply(Matrix_Sqrt_Apply(v)) ~= Matrix_Apply(v)
            % because the preconditioned system is a factor, i.e., A=S*S^T, but not a square root
        end

        function [x12, relres] = Krylov_Sqrt(this, A, b, maxiter, tol)
            % This function computes A^{1/2}b using Lanczos approach
            % described in Algorithm 2.2 of
            %   "Quantifying Uncertainties in Bayesian Linear Inverse Problems using
            %   Krylov Subspace Methods" - Saibaba, Chung, and Petroske, 2018
            %
            % Inputs:
            %   A (n x n) - func type
            %   b (n x 1) - right hand side
            %   maxiter   - maximum number of Lanczos iterations
            %   tol - tolerance for stopping

            n = size(b, 1);
            nrmb = norm(b);

            % Initialize Lanczos quantities
            V = zeros(n, maxiter);
            T = zeros(maxiter + 1, maxiter + 1);

            % First step
            vj = b / nrmb;
            vjm1 = b * 0;
            beta = 0;
            ykp = 0;
            relres = zeros(maxiter, 1);

            for j = 1:maxiter
                V(:, j) = vj;
                wj = A(vj);
                alpha = wj' * vj;
                wj = wj - alpha * vj - beta * vjm1;
                beta = norm(wj);

                % Set vectors for new iterations
                vjm1 = vj;
                vj =    wj / beta;

                % Reorthogonalize vj (CGS2)
                vj = vj - V(:, 1:j) * (V(:, 1:j)' * vj);
                vj = vj / norm(vj);
                vj = vj - V(:, 1:j) * (V(:, 1:j)' * vj);
                vj = vj / norm(vj);

                % Set the tridiagonal matrix
                T(j, j) = alpha;
                T(j + 1, j) = beta;
                T(j, j + 1) = beta;

                % Compute partial Lanczos solution
                Tk = T(1:j, 1:j);
                Vk = V(:, 1:j);
                e1 = zeros(j, 1);
                e1(1) = nrmb;
                T12 = sqrtm(Tk);
                yk = T12 * e1;

                %   relres(j) = norm(A*(Vk*(Tk\e1))))-b); % Lanczos residual
                % Check differences b/w successive iterations
                relres(j) = norm([ykp; 0] - yk) / norm(yk);
                if  relres(j) < tol
                    relres = relres(1:j);
                    break
                else
                    ykp = yk;
                end
            end
            x12 = Vk * yk;
        end

    end

end
