function correlation_length = computeCorrelationLength_1D(x,u)

xl = linspace(min(x), max(x), 100)';
f = interp1(x,u,xl);

% Compute the mean and mean square of the function
f_mean = mean(f(:));
f_mean_sq = mean(f(:).^2);

r_samples = 30;
max_r = max(x)-min(x);
r_range = linspace(0,max_r,r_samples)';

num_samples = 100;
line_samples = linspace(-1,1,num_samples)';

C = zeros(r_samples,1);
index = r_samples;
for r = 1:r_samples
    mean_r = zeros(num_samples,1);
    for k = 1:num_samples
        xk = xl + line_samples(k)*r_range(r);
        xk(xk>max(xl)) = NaN;
        xk(xk<min(xl)) = NaN;
        fk = interp1(x,u,xk);
        mean_r(k) = mean((f-f_mean).*(fk-f_mean),'omitnan');
    end
    C(r) = mean(mean_r)/(f_mean_sq-f_mean^2);
    if C(r) < C(1)/exp(1)
        index = r;
        break
    end
end

correlation_length = r_range(index);
end