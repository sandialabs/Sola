function correlation_length = computeCorrelationLength_2D(x,y,u)

xl = linspace(min(x), max(x), 100)';
yl = linspace(min(y), max(y), 100)';
[X, Y] = meshgrid(xl, yl);
F = scatteredInterpolant(x, y, u);
f = F(X, Y);

% Compute the mean and mean square of the function
f_mean = mean(f(:));
f_mean_sq = mean(f(:).^2);

r_samples = 30;
max_r = norm([min(x);min(y)]-[max(x);max(y)]);
r_range = linspace(0,max_r,r_samples)';

num_samples = 100;
circle_samples = this.sampleUnitCircle(num_samples);

C = zeros(r_samples,1);
index = r_samples;
for r = 1:r_samples
    mean_r = zeros(num_samples,1);
    for k = 1:num_samples
        Xk = X + circle_samples(k,1)*r_range(r);
        Yk = Y + circle_samples(k,2)*r_range(r);
        Xk(Xk>max(xl)) = NaN;
        Xk(Xk<min(xl)) = NaN;
        Yk(Yk>max(yl)) = NaN;
        Yk(Yk<min(yl)) = NaN;
        fk = F(Xk,Yk);
        mean_r(k) = mean((f(:)-f_mean).*(fk(:)-f_mean),'omitnan');
    end
    C(r) = mean(mean_r)/(f_mean_sq-f_mean^2);
    if C(r) < C(1)/exp(1)
        index = r;
        break
    end
end

correlation_length = r_range(index);
end