classdef Example_2_Objective < Dynamic_Objective
    % Solve the optimiation problem
    % min_{z} J(S(z),z) = int_0^T g(S(z)(t))dt + R(z)
    % where S(z) solves the odinary differential equation
    % dy/dt = [y_1 , y_2]
    % y(0) = [z_1 , z_2]
    % g(y) = (y_1-exp(t))^2 + (y_2-exp(t))^2
    % R(z) = (z_1-1)^2 + (z_2-1)^2

    properties

    end

    methods (Access = public)

        function [val, grad_y] = Time_Instance_Objective(this, y, t)
            val = (y(1) - exp(t))^2 + (y(2) - exp(t))^2;
            grad_y = 2 * (y - exp(t));
        end

        function [val, grad_z] = Regularization_Objective(this, z)
            val = (z(1) - 1)^2 + (z(2) - 1)^2;
            grad_z = 2 * (z - 1);
        end

        function [Mv] = Time_Instance_Objective_yy_Apply(this, v, y, t)
            A = 2 * eye(2);
            Mv = A * v;
        end

        function [Mv] = Regularization_Objective_zz_Apply(this, v, z)
            A = 2 * eye(2);
            Mv = A * v;
        end

    end

    methods (Access = public)

        function this = Example_2_Objective(m, n, T, N)
            this = this@Dynamic_Objective(m, n, T, N);
        end

    end

end
