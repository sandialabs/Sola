%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Tutorial_2B_Constraint < Dynamic_Constraint
    % Constraint for an orbital dynamics problem.
    %
    % dy/dt = (y2, y1 y4^2 - k/y1^2 + q2, y4, (q1 - 2 y2 y4) / y1)
    % y(0) = z
    %
    % Here y = (r, r', w, w') represents radius, change in radius, angular
    % momentum, and change in angular momentum of an object. The constants
    % r and w are the target radius and angular momentum. The controls
    % (q1(t), q2(t)) are thrust from the satellite in the radial and
    % tangential directions, respectively.

    properties
        k            % Constant of proportionality.
        n_q          % Number of controls at each fixed time.
        y0           % Initial condition.
    end

    methods (Access = public)

        function this = Tutorial_2B_Constraint(T, n_t, r0, w0, k)
            % Constructor.
            %
            % Parameters
            % ----------
            % T
            %   Final simulation time.
            % n_t
            %   Number of time steps.
            % r0
            %   Initial radius of the satellite.
            % w0
            %   Initial angular velocity of the satellite.
            n_q = 2;
            n_z = n_q * (n_t - 1);
            this = this@Dynamic_Constraint(4, n_z, T, n_t);
            this.y0 = [r0; 0; 0; w0];
            this.k = k;
            this.n_q = n_q;
        end

        function [f, f_y, f_z] = f(this, y, z, t)
            I = this.Input_Indices(t);
            q = z(I);

            f = [y(2)
                 y(1) * y(4)^2 - (this.k / y(1)^2) + q(1)
                 y(4)
                 (q(2) - 2 * y(2) * y(4)) / y(1)];

            f_y = [0, 1, 0, 0
                   y(4)^2 + 2 * (this.k / y(1)^3), 0, 0, 2 * y(1) * y(4)
                   0, 0, 0, 1
                   (2 * y(2) * y(4) - q(2)) / y(1)^2, (-2 * y(4)) / y(1), 0, (-2 * y(2)) / y(1)];

            f_z = zeros(this.n_y, this.n_z);
            f_z(:, I) = [0, 0
                         1, 0
                         0, 0
                         0, 1 / y(1)];
        end

        function [h, h_z] = h(this, ~)
            h = this.y0;
            h_z = zeros(this.n_y, this.n_z);
        end

        function [y_out] = f_yy_Apply(this, y_in, y, z, t, lambda)
            I = this.Input_Indices(t);
            q = z(I);

            ijk211 = -6 * this.k * lambda(2) * y_in(1, :) / y(1)^4;
            ijk214 = 2 * lambda(2) * y(4) * y_in(4, :);
            ijk241 = 2 * lambda(2) * y(4) * y_in(1, :);
            ijk244 = 2 * lambda(2) * y(1) * y_in(4, :);
            ijk411 = 2 * lambda(4) * y_in(1, :) * (q(2) - 2 * y(2) * y(4)) / y(1)^3;
            ijk412 = 2 * lambda(4) * y(4) * y_in(2, :) / y(1)^2;
            ijk414 = 2 * lambda(4) * y(2) * y_in(4, :) / y(1)^2;
            ijk421 = 2 * lambda(4) * y(4) * y_in(1, :) / y(1)^2;
            ijk441 = 2 * lambda(4) * y(2) * y_in(1, :) / y(1)^2;

            y_out = [ijk211 + ijk214 + ijk411 + ijk412 + ijk414
                     ijk421
                     zeros(1, size(y_in, 2))
                     ijk241 + ijk244 + ijk441];
        end

        function [y_out] = f_yz_Apply(this, z_in, y, ~, ~, lambda)
            nvecs = size(z_in, 2);
            y_out = zeros(this.n_y, nvecs);
            mask = 2:2:this.n_z;
            y_out(1, :) = sum(-lambda(4) * z_in(mask, :) / y(1)^2);
        end

        function [z_out] = f_zy_Apply(this, y_in, y, ~, ~, lambda)
            z_out = zeros(this.n_z, 1);
            mask = 2:2:this.n_z;
            z_out(mask) = -lambda(4) * y_in(1) / y(1)^2;
        end

        function [z_out] = f_zz_Apply(this, z_in, ~, ~, ~, ~)
            z_out = zeros(this.n_z, size(z_in, 2));
        end

        function [z_out] = h_zz_Apply(this, z_in, ~, ~)
            z_out = zeros(this.n_z, size(z_in, 2));
        end

    end

    methods (Access = private)

        function [mask] = Input_Indices(this, t)
            % Get the indices of the control at time t, i.e.,
            % I = Input_Indices(t) --> q(t) = z(I).

            [~, t_index] = min(abs(t - this.t_mesh));
            if t_index == 1
                error('no control at initial time!');
            end
            idx = t_index - 1;
            mask = (this.n_q * (idx - 1) + 1):(this.n_q * idx);
        end

    end
end
