classdef Example_2_Objective < Dynamic_Objective
    % Solve the optimization problem
    % min_{z} J(S(z),z) = int_0^T g(S(z)(t))dt + R(z)
    % where S(z) solves the ordinary differential equation
    % dy/dt = [y_1 , y_2]
    % y(0) = [z_1 , z_2]
    % g(y) = (y_1-exp(t))^2 + (y_2-exp(t))^2
    % R(z) = (z_1-1)^2 + (z_2-1)^2

    properties

    end

    methods (Access = public)

        function [val, grad_y] = g(this, y, t)
            val = (y(1) - exp(t))^2 + (y(2) - exp(t))^2;
            grad_y = 2 * (y - exp(t));
        end

        function [val, grad_z] = R(this, z)
            val = (z(1) - 1)^2 + (z(2) - 1)^2;
            grad_z = 2 * (z - 1);
        end

        function [Mv] = g_yy_Apply(this, v, y, t)
            Mv = 2 * eye(2) * v;
        end

        function [Mv] = R_zz_Apply(this, v, z)
            Mv = 2 * eye(2) * v;
        end

    end

    methods (Access = public)

        function this = Example_2_Objective(n_y, n_z, T, n_t)
            this = this@Dynamic_Objective(n_y, n_z, T, n_t);
        end

    end

end
