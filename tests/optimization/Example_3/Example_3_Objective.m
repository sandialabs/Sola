classdef Example_3_Objective < Dynamic_Objective

    % Solve the optimiation problem
    % min_{z} J(S(z),z) = int_0^T g(S(z)(t))dt + R(z)
    % where S(z) solves the odinary differential equation
    % dy/dt = [y3/y2, 2*y1, 3*y1*y2]
    % y(0) = [z_1^2 , z_2^3, z_3^4]
    % g(y) = (y_1-exp(t))^2 + (y_2-exp(2*t))^2 + (y_3-exp(3*t))^2
    % R(z) = (z_1-1)^2 + (z_2-1)^2 + (z_3-1)^2

    methods (Access = public)

        function [val, grad_y] = g(this, y, t)
            val = (y(1) - exp(t))^2 + (y(2) - exp(2 * t))^2 + (y(3) - exp(3 * t))^2;
            grad_y = zeros(3, 1);
            grad_y(1) = 2 * (y(1) - exp(t));
            grad_y(2) = 2 * (y(2) - exp(2 * t));
            grad_y(3) = 2 * (y(3) - exp(3 * t));
        end

        function [val, grad_z] = R(this, z)
            val = (z(1) - 1)^2 + (z(2) - 1)^2 + (z(3) - 1)^2;
            grad_z = 2 * (z - 1);
        end

        function [Mv] = g_yy_Apply(this, v, y, t)
            Mv = 2 * eye(3) * v;
        end

        function [Mv] = R_zz_Apply(this, v, z)
            Mv = 2 * eye(3) * v;
        end

    end

    methods (Access = public)

        function this = Example_3_Objective(n_y, n_z, T, n_t)
            this = this@Dynamic_Objective(n_y, n_z, T, n_t);
        end

    end

end
