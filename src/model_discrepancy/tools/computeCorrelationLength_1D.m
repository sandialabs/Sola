function correlation_length = computeCorrelationLength_1D(x,u)

num_samples = 100;
xl = linspace(min(x), max(x), num_samples)';
f = interp1(x,u,xl);

% Compute the mean and variance of the function
f_mean = mean(f(:));
f_var = (1/(num_samples-1))*sum((f-f_mean).^2);

r_samples = 30;
max_r = max(x)-min(x);
r_range = linspace(0,max_r,r_samples+1)';
r_range = r_range(2:end);
line_samples = linspace(-1,1,num_samples)';

C = zeros(r_samples,1);
index = nan;
for r = 1:r_samples
    cov_r = zeros(num_samples,1);
    for k = 1:num_samples
        xk = xl + line_samples(k)*r_range(r);
        xk(xk>max(xl)) = NaN;
        xk(xk<min(xl)) = NaN;
        fk = interp1(x,u,xk);
        Nk = sum(0*fk+1,'omitnan');
        cov_r(k) = (1/(Nk-1)) * sum((f-f_mean).*(fk-f_mean),'omitnan');
    end
    C(r) = mean(cov_r)/f_var;
    if C(r) < .1 
        index = r;
        break
    end
end

if isnan(index)
    correlation_length = nan;
else
    correlation_length = r_range(index);
end
end