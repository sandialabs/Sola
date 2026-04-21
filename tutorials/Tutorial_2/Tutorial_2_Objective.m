%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Tutorial_2_Objective < Dynamic_Objective
    % Objective function for an orbital dynamics problem.
    %
    % min \int_{0}^{T} ||y - (r, 0, wt, t)^T)||^2 dt.
    %
    % Here y = (r, r', w, w') represents radius, change in radius, angular
    % momentum, and change in angular momentum of an object. The constants
    % r and w are the target radius and angular velocity.

    properties
        radius       % Target orbital radius.
        velocity     % Target angular velocity.
    end

    methods

        function this = Tutorial_2_Objective(T, n_t, radius, velocity)
            this = this@Dynamic_Objective(4, 2, T, n_t);
            this.radius = radius;
            this.velocity = velocity;
        end

        function [val, grad_y] = g(this, y, t)
            w = this.velocity;
            r = this.radius;
            a = [r; 0; w * t; w];
            y_minus_a = y - a;

            val = sum(y_minus_a.^2);
            grad_y = 2 * y_minus_a;
        end

        function [val, grad_z] = R(this, ~)
            val = 0;
            grad_z = zeros(this.n_z, 1);
        end

        function [y_out] = g_yy_Apply(~, y_in, ~, ~)
            y_out = 2 * y_in;
        end

        function [z_out] = R_zz_Apply(this, z_in, ~)
            z_out = zeros(this.n_z, size(z_in, 2));
        end

        %% Plotting functions.
        function Plot(this, u)
            % Plot a trajectory two ways: x(t) and y(t) in time, and
            % (x(t),y(t)) as a 2D trajectory.
            %
            % Parameters
            % ----------
            % u
            %   Solution array, either an :math:`n_y \times n_t` matrix or
            %   an :math:`(n_y n_t, 1)` vector.

            t = this.t_mesh;
            x_goal = this.radius * cos(this.velocity .* t);
            y_goal = this.radius * sin(this.velocity .* t);

            y = reshape(u, this.n_y, this.n_t);
            r = y(1, :);
            theta = y(3, :);
            x1 = r .* cos(theta);
            x2 = r .* sin(theta);
            lim = 1.2 * max(r);

            % Plot the state trajectory coordinates in time.
            fig = figure;
            fig.Position(3:4) = [830, 300];
            subplot(1, 2, 1);
            plot(t, x1, '-', 'LineWidth', 2);
            hold on;
            plot(t, x2, '-', 'LineWidth', 2);
            plot(t, -ones(this.n_t, 1), 'k-', 'LineWidth', 0.1);
            plot(t, ones(this.n_t, 1), 'k-', 'LineWidth', 0.1);
            xlim([0, t(end)]);
            ylim([-lim, lim]);
            xlabel('$$t$$', 'Interpreter', 'latex');
            title('Position coordinates');
            legend({'$$x_1(t)$$', '$$x_2(t)$$', '', ''}, ...
                   'Location', 'southeast', 'Interpreter', 'latex');

            % Plot the state trajectory in two-dimensional space.
            subplot(1, 2, 2);
            plot(x_goal, y_goal, 'k--', 'LineWidth', 2);
            hold on;
            plot(x1, x2, '-', 'LineWidth', 1);
            plot(x1(1), x2(1), '.', 'MarkerSize', 16);
            xlim([-lim, lim]);
            ylim([-lim, lim]);
            xlabel('$$x_1(t)$$', 'Interpreter', 'latex');
            ylabel('$$x_2(t)$$', 'Interpreter', 'latex');
            title('Position');
            axis('equal');
            legend({'target trajectory', 'realized trajectory', 'initial position'}, ...
                   'Location', 'southeast', 'Interpreter', 'latex');
        end

    end
end
