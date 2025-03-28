function correlation_length = computeCorrelationLength_2D(x, y, u, initial_guess)

    N = 100;
    xl = linspace(min(x), max(x), N)';
    yl = linspace(min(y), max(y), N)';
    [X, Y] = meshgrid(xl, yl);
    F = scatteredInterpolant(x, y, u);
    f = F(X, Y);

    % Compute the mean and variance of the function
    f_mean = mean(f(:));
    f_var = (1 / (N^2 - 1)) * sum((f(:) - f_mean).^2);

    r_samples = 60;
    max_r = norm([min(x); min(y)] - [max(x); max(y)]);
    r_range = linspace(0, max_r, r_samples + 1)';
    r_range = r_range(2:end);

    [~, r_init] = min(abs(r_range - initial_guess));
    index = r_init;

    num_samples = 20;
    circle_samples = sampleUnitCircle(num_samples);

    cov_r = zeros(num_samples, 1);
    for k = 1:num_samples
        Xk = X + circle_samples(k, 1) * r_range(r_init);
        Yk = Y + circle_samples(k, 2) * r_range(r_init);
        Xk(Xk > max(xl)) = NaN;
        Xk(Xk < min(xl)) = NaN;
        Yk(Yk > max(yl)) = NaN;
        Yk(Yk < min(yl)) = NaN;
        fk = F(Xk, Yk);
        Nk = sum(0 * fk(:) + 1, 'omitnan');
        cov_r(k) = (1 / (Nk - 1)) * sum((f(:) - f_mean) .* (fk(:) - f_mean), 'omitnan');
    end
    C = mean(cov_r) / f_var;

    if C > .1
        for r = (r_init + 1):r_samples
            cov_r = zeros(num_samples, 1);
            for k = 1:num_samples
                Xk = X + circle_samples(k, 1) * r_range(r);
                Yk = Y + circle_samples(k, 2) * r_range(r);
                Xk(Xk > max(xl)) = NaN;
                Xk(Xk < min(xl)) = NaN;
                Yk(Yk > max(yl)) = NaN;
                Yk(Yk < min(yl)) = NaN;
                fk = F(Xk, Yk);
                Nk = sum(0 * fk(:) + 1, 'omitnan');
                cov_r(k) = (1 / (Nk - 1)) * sum((f(:) - f_mean) .* (fk(:) - f_mean), 'omitnan');
            end
            C = mean(cov_r) / f_var;
            if C < .1
                index = r;
                break
            end
        end
    else
        for r = (r_init - 1):-1:1
            cov_r = zeros(num_samples, 1);
            for k = 1:num_samples
                Xk = X + circle_samples(k, 1) * r_range(r);
                Yk = Y + circle_samples(k, 2) * r_range(r);
                Xk(Xk > max(xl)) = NaN;
                Xk(Xk < min(xl)) = NaN;
                Yk(Yk > max(yl)) = NaN;
                Yk(Yk < min(yl)) = NaN;
                fk = F(Xk, Yk);
                Nk = sum(0 * fk(:) + 1, 'omitnan');
                cov_r(k) = (1 / (Nk - 1)) * sum((f(:) - f_mean) .* (fk(:) - f_mean), 'omitnan');
            end
            C = mean(cov_r) / f_var;
            if C < .1
                index = r;
                break
            end
        end
    end

    correlation_length = r_range(index);
end
