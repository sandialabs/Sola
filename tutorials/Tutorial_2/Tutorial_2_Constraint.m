%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Tutorial_2_Constraint < Dynamic_Constraint
    % Constraint for an orbital dynamics problem.
    %
    % dy/dt = (y2, y1 y4^2 - k/y1^2, y4, -2 y2 y4 / y1)\trp
    % y(0) = (z_1, 0, 0, z_2)\trp
    %
    % Here y = (r, r', theta, theta') represents radius, change in radius,
    % angular momentum, and change in angular momentum of an object.
    % The constant k accounts for the mass of the object being orbited and
    % a constant of proportionality.

    properties
        k            % Constant of proportionality.
    end

    methods (Access = public)

        function this = Tutorial_2_Constraint(T, n_t, k)
            this = this@Dynamic_Constraint(4, 2, T, n_t);
            this.k = k;
        end

        function [f, f_y, f_z] = f(this, y, ~, ~)
            f = [y(2)
                 y(1) * y(4)^2 - (this.k / y(1)^2)
                 y(4)
                 -2 * y(2) * y(4) / y(1)];
            f_y = [0, 1, 0, 0
                   y(4)^2 + 2 * (this.k / y(1)^3), 0, 0, 2 * y(1) * y(4)
                   0, 0, 0, 1
                   2 * y(2) * y(4) / y(1)^2, (-2 * y(4)) / y(1), 0, (-2 * y(2)) / y(1)];

            f_z = zeros(this.n_y, this.n_z);
        end

        function [h, h_z] = h(~, z)
            h = [z(1); 0; 0; z(2)];
            h_z = [1 0
                   0 0
                   0 0
                   0 1];
        end

        function [y_out] = f_yy_Apply(this, y_in, y, ~, ~, lambda)
            ijk211 = -6 * this.k * lambda(2) * y_in(1, :) / y(1)^4;
            ijk214 = 2 * lambda(2) * y(4) * y_in(4, :);
            ijk241 = 2 * lambda(2) * y(4) * y_in(1, :);
            ijk244 = 2 * lambda(2) * y(1) * y_in(4, :);
            ijk411 = -4 * lambda(4) * y(2) * y(4) * y_in(1, :) / y(1)^3;
            ijk412 = 2 * lambda(4) * y(4) * y_in(2, :) / y(1)^2;
            ijk414 = 2 * lambda(4) * y(2) * y_in(4, :) / y(1)^2;
            ijk424 = -2 * lambda(4) * y_in(4, :) / y(1);
            ijk421 = 2 * lambda(4) * y(4) * y_in(1, :) / y(1)^2;
            ijk441 = 2 * lambda(4) * y(2) * y_in(1, :) / y(1)^2;
            ijk442 = -2 * lambda(4) * y_in(2, :) / y(1);
            y_out = [ijk211 + ijk214 + ijk411 + ijk412 + ijk414
                     ijk421 + ijk424
                     zeros(1, size(y_in, 2))
                     ijk241 + ijk244 + ijk441 + ijk442];
        end

        function [y_out] = f_yz_Apply(this, z_in, ~, ~, ~, ~)
            y_out = zeros(this.n_y, size(z_in, 2));
        end

        function [z_out] = f_zy_Apply(this, y_in, ~, ~, ~, ~)
            z_out = zeros(this.n_z, size(y_in, 2));
        end

        function [z_out] = f_zz_Apply(this, z_in, ~, ~, ~, ~)
            z_out = zeros(this.n_z, size(z_in, 2));
        end

        function [z_out] = h_zz_Apply(this, z_in, ~, ~)
            z_out = zeros(this.n_z, size(z_in, 2));
        end

    end
end
