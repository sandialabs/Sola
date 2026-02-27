function [x12, relres] = krylov_sqrt(A, b, maxiter, tol)
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
