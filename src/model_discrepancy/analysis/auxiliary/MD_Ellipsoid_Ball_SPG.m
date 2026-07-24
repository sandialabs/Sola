classdef MD_Ellipsoid_Ball_SPG
    %MD_ELLIPSOID_BALL_SPG
    %
    % Utility class for minimizing a smooth function over a product of
    % weighted ellipsoid balls using a spectral projected gradient method.
    %
    % The feasible set is
    %
    %   ||B(:, j) - beta_bar||_A <= radius,     j = 1, ..., p,
    %
    % equivalently,
    %
    %   (B(:, j) - beta_bar)' A (B(:, j) - beta_bar) <= radius^2,
    %   j = 1, ..., p,
    %
    % where A = Q * diag(d) * Q'.

    methods (Static)

        function projection_data = Prepare_Projection(beta_bar, p, constr_radius, Q, d)

            beta_bar = beta_bar(:);
            r = length(beta_bar);

            if p < 1 || p ~= round(p)
                error('p must be a positive integer.');
            end

            if constr_radius < 0
                error('The ellipsoid constraint radius must be nonnegative.');
            end

            if size(Q, 1) ~= r || size(Q, 2) ~= r
                error('Q must be an r-by-r matrix, where r = length(beta_bar).');
            end

            d = d(:);

            if length(d) ~= r
                error('Length of d must equal length(beta_bar).');
            end

            if any(d < -1.e-14 * max(1, max(abs(d))))
                error('Ellipsoid eigenvalues must be nonnegative.');
            end

            % Remove tiny negative eigenvalues caused by roundoff.
            d(d < 0) = 0;

            projection_data = struct;
            projection_data.r = r;
            projection_data.p = p;
            projection_data.radius = constr_radius;

            projection_data.beta_bar = beta_bar;
            projection_data.Q = Q;
            projection_data.d = d;
            projection_data.QT_beta_bar = Q.' * beta_bar;

        end

        function [x, info] = Minimize(fun, x0, projection_data, opts)

            if nargin < 4
                opts = struct;
            end

            max_iter = MD_Ellipsoid_Ball_SPG.Get_Option(opts, 'max_iter', 1000);
            pg_tol = MD_Ellipsoid_Ball_SPG.Get_Option(opts, 'pg_tol', 1.e-8);
            armijo_c = MD_Ellipsoid_Ball_SPG.Get_Option(opts, 'armijo_c', 1.e-4);
            backtrack_factor = MD_Ellipsoid_Ball_SPG.Get_Option(opts, 'backtrack_factor', 0.5);
            max_backtracks = MD_Ellipsoid_Ball_SPG.Get_Option(opts, 'max_backtracks', 60);
            nonmonotone_window = MD_Ellipsoid_Ball_SPG.Get_Option(opts, 'nonmonotone_window', 10);
            verbosity = MD_Ellipsoid_Ball_SPG.Get_Option(opts, 'verbosity', false);

            x = MD_Ellipsoid_Ball_SPG.Project(x0, projection_data);

            [f, g] = fun(x);

            alpha = 1;

            f_hist = nan(max_iter + 1, 1);
            f_hist(1) = f;
            hist_len = 1;

            converged = false;
            pg_norm = NaN;
            iter = 0;

            n_f_eval = 0;
            n_fg_eval = 1;

            for iter = 1:max_iter

                x_pg = MD_Ellipsoid_Ball_SPG.Project(x - g, projection_data);
                pg = x - x_pg;
                pg_norm = norm(pg);

                if pg_norm <= pg_tol * max(1, norm(x))
                    converged = true;
                    break;
                end

                % Reuse x_pg when alpha == 1.
                if alpha == 1
                    x_proj = x_pg;
                else
                    x_proj = MD_Ellipsoid_Ball_SPG.Project(x - alpha * g, projection_data);
                end

                dir = x_proj - x;
                gtd = g.' * dir;

                if gtd >= 0
                    alpha = 1;
                    x_proj = x_pg;
                    dir = x_proj - x;
                    gtd = g.' * dir;

                    if gtd >= 0
                        break;
                    end
                end

                hist_start = max(1, hist_len - nonmonotone_window + 1);
                f_ref = max(f_hist(hist_start:hist_len));

                accepted = false;
                step = 1;

                for bt = 1:max_backtracks

                    x_trial = x + step * dir;
                    [f_trial, g_trial] = fun(x_trial);
                    n_fg_eval = n_fg_eval + 1;

                    if f_trial <= f_ref + armijo_c * step * gtd
                        accepted = true;
                        break;
                    end

                    step = backtrack_factor * step;
                end

                if ~accepted
                    x_trial = x + step * dir;
                    [f_trial, g_trial] = fun(x_trial);
                    n_fg_eval = n_fg_eval + 1;
                end

                s = x_trial - x;
                y = g_trial - g;

                x = x_trial;
                f = f_trial;
                g = g_trial;

                hist_len = hist_len + 1;
                f_hist(hist_len) = f;

                % Barzilai-Borwein spectral step.
                sy = s.' * y;

                if sy > eps * norm(s) * max(1, norm(y))
                    alpha = (s.' * s) / sy;
                    alpha = min(max(alpha, 1.e-12), 1.e12);
                else
                    alpha = 1;
                end

            end

            info = struct;
            info.iter = iter;
            info.final_objective = f;
            info.projected_gradient_norm = pg_norm;
            info.converged = converged;
            info.f_hist = f_hist(1:hist_len);
            info.n_f_eval = n_f_eval;
            info.n_fg_eval = n_fg_eval;

            if verbosity
                fprintf('SPG completed in %d iterations. Final objective: %.6e\n', iter, f);
                fprintf('Projected gradient norm: %.6e\n', pg_norm);
                fprintf('Objective-only evaluations: %d\n', n_f_eval);
                fprintf('Objective/gradient evaluations: %d\n', n_fg_eval);
            end

        end

        function x = Project(y, projection_data)
            %PROJECT Euclidean projection onto product of ellipsoid balls.
            %
            % Projects y onto the feasible set
            %
            %   ||B(:, j) - beta_bar||_A <= radius,     j = 1, ..., p,
            %
            % where A = Q * diag(d) * Q'.

            r = projection_data.r;
            p = projection_data.p;
            Q = projection_data.Q;
            d = projection_data.d;
            radius_sq = projection_data.radius^2;

            y = y(:);

            if mod(length(y), r) ~= 0
                error('Input vector length must be divisible by the reduced dimension r.');
            end

            if length(y) ~= r * p
                error('Input vector length is inconsistent with the projection data.');
            end

            Y = reshape(y, r, p);
            X = zeros(r, p);
            Y_hat = Q.' * Y - projection_data.QT_beta_bar;

            for j = 1:p

                y_hat_j = Y_hat(:, j);

                current_radius_sq = sum(d .* y_hat_j.^2);

                if current_radius_sq <= radius_sq * (1 + 1.e-12)
                    % Already feasible.
                    x_hat_j = y_hat_j;
                else
                    % Euclidean projection onto
                    %
                    %   sum_i d_i * x_i^2 <= radius^2.
                    %
                    % The projected point has the form
                    %
                    %   x_i = y_i / (1 + 2 * lambda * d_i),
                    %
                    % where lambda >= 0 is chosen so that the constraint is active.
                phi = @(lambda) ...
                        sum(d .* y_hat_j.^2 ./ (1 + 2 * lambda * d).^2) ...
                        - radius_sq;

                    lambda_lo = 0;
                    lambda_hi = 1;

                    while phi(lambda_hi) > 0
                        lambda_hi = 2 * lambda_hi;
                end

                    for k = 1:60
                    lambda_mid = 0.5 * (lambda_lo + lambda_hi);

                        if phi(lambda_mid) > 0
                            lambda_lo = lambda_mid;
                        else
                            lambda_hi = lambda_mid;
                end
                    end

                    lambda = lambda_hi;

                    x_hat_j = y_hat_j ./ (1 + 2 * lambda * d);
            end

                X(:, j) = projection_data.beta_bar + Q * x_hat_j;

            end

            x = X(:);
        end

    end
    
    methods (Static, Access = private)

        function value = Get_Option(opts, name, default_value)

            if isfield(opts, name) && ~isempty(opts.(name))
                value = opts.(name);
            else
                value = default_value;
            end

        end

    end
end