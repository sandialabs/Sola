classdef Tutorial_2_Constraint_AD < Dynamic_Constraint_AD
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
        k
    end

    methods (Access = public)

        function this = Tutorial_2_Constraint_AD(T, n_t, k)
            this = this@Dynamic_Constraint_AD(4, 2, T, n_t);
            this.k = k;
        end

        function [f] = f_AD(this, y, ~, ~)
            f = [y(2)
                 y(1) * y(4)^2 - (this.k / y(1)^2)
                 y(4)
                 -2 * y(2) * y(4) / y(1)];
        end

        function [h] = h_AD(~, z)
            h = [z(1); 0; 0; z(2)];
        end

    end
end
