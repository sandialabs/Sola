classdef Example_4_Objective < Dynamic_Objective

    % Solve the optimiation problem
    % min_{z} J(S(z),z) = int_0^T g(S(z)(t))dt + R(z)
    % where S(z) solves the odinary differential equation
    % dy/dt = [z_1^2*y3/y2, 2*z_2^2*y1, 3*z_3^2*y1*y2]
    % y(0) = [1, 1, 1]
    % g(y) = (y_1-exp(t))^2 + (y_2-exp(2*t))^2 + (y_3-exp(3*t))^2
    % R(z) = (z_1-1)^2 + (z_2-1)^2 + (z_3-1)^2

    methods (Access = public)

        function [val, grad_y] = Time_Instance_Objective(this, y, t)
            val = (y(1) - exp(t))^2 + (y(2) - exp(2 * t))^2 + (y(3) - exp(3 * t))^2;
            grad_y = zeros(3, 1);
            grad_y(1) = 2 * (y(1) - exp(t));
            grad_y(2) = 2 * (y(2) - exp(2 * t));
            grad_y(3) = 2 * (y(3) - exp(3 * t));
        end

        function [val, grad_z] = Regularization_Objective(this, z)
            val = (z(1) - 1)^2 + (z(2) - 1)^2 + (z(3) - 1)^2;
            grad_z = 2 * (z - 1);
        end

        function [Mv] = Time_Instance_Objective_yy_Apply(this, v, y, t)
            A = 2 * eye(3);
            Mv = A * v;
        end

        function [Mv] = Regularization_Objective_zz_Apply(this, v, z)
            A = 2 * eye(3);
            Mv = A * v;
        end

    end

end
