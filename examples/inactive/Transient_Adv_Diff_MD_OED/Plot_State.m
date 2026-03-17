function [] = Plot_State(x, u)
    n_y = length(x);
    n_t = length(u) / n_y;
    ymin = min(u) - .05 * abs(min(u));
    ymax = max(u) + .05 * abs(max(u));
    u_reshape = reshape(u, n_y, n_t);
    figure;
    hold on;
    for j = 1:n_t
        t = (j - 1) / (n_t - 1);
        plot(x, u_reshape(:, j), 'LineWidth', 3, 'Color', (1 - .8 * t) * ones(3, 1));
        ylim([ymin, ymax]);
        title(['Time = ', num2str(t)]);
        pause(.05);
    end
end
