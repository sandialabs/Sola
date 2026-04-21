%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function correlation_length = Compute_Correlation_Length_1D(x, u, initial_guess)

    N = size(u, 2);
    num_samples = 200;
    xl = linspace(min(x), max(x), num_samples)';
    f = interp1(x, u, xl);

    % Compute the mean and variance of the function
    f_mean = mean(f(:));
    f_var = (1 / (N * num_samples - 1)) * sum((f(:) - f_mean).^2);

    [~, r_init] = min(abs(xl - initial_guess));
    index = r_init;
    cov_r = (1 / (N * (num_samples - r_init) - 1)) * sum(sum((f(1:(end - r_init), :) - f_mean) .* (f((1 + r_init):end, :) - f_mean)));
    C = cov_r / f_var;
    if C > .1
        for r = (r_init + 1):(num_samples - 1)
            cov_r = (1 / (N * (num_samples - r) - 1)) * sum(sum((f(1:(end - r), :) - f_mean) .* (f((1 + r):end, :) - f_mean)));
            C = cov_r / f_var;
            if C < .1
                index = r;
                break
            end
        end
    else
        for r = (r_init - 1):-1:1
            cov_r = (1 / (N * (num_samples - r) - 1)) * sum(sum((f(1:(end - r), :) - f_mean) .* (f((1 + r):end, :) - f_mean)));
            C = cov_r / f_var;
            if C > .1
                index = r - 1;
                break
            end
        end
    end

    correlation_length = xl(index + 1) - xl(1);

end
