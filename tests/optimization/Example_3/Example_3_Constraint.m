classdef Example_3_Constraint < Dynamic_Constraint

    % Solve the optimization problem
    % min_{z} J(S(z),z) = int_0^T g(S(z)(t))dt + R(z)
    % where S(z) solves the ordinary differential equation
    % dy/dt = [y3/y2, 2*y1, 3*y1*y2]
    % y(0) = [z_1^2 , z_2^3, z_3^4]
    % g(y) = (y_1-exp(t))^2 + (y_2-exp(2*t))^2 + (y_3-exp(3*t))^2
    % R(z) = (z_1-1)^2 + (z_2-1)^2 + (z_3-1)^2

    methods (Access = public)

        function [f, f_y, f_z] = Time_Instance_RHS(this, y, z, t)
            f = [y(3) / y(2); 2 * y(1); 3 * y(1) * y(2)];
            f_y = [0, -y(3) / y(2)^2, 1 / y(2); 2, 0, 0; 3 * y(2), 3 * y(1), 0];
            f_z = zeros(3, 3);
        end

        function [h, h_z] = Initial_Condition(this, z)
            h = [z(1)^2; z(2)^3; z(3)^4];
            h_z = [2 * z(1), 0, 0; 0, 3 * z(2).^2, 0; 0, 0, 4 * z(3)^3];
        end

        function [Mv] = Time_Instance_RHS_yy_Apply(this, v, y, z, t, lambda)
            A = [0, 3 * lambda(3), 0; 3 * lambda(3), 2 * lambda(1) * y(3) / y(2)^3, -lambda(1) / y(2)^2; 0, -lambda(1) / y(2)^2, 0];
            Mv = A * v;
        end

        function [Mv] = Time_Instance_RHS_yz_Apply(this, v, y, z, t, lambda)
            Mv = zeros(this.n_y, size(v, 2));
        end

        function [Mv] = Time_Instance_RHS_zy_Apply(this, v, y, z, t, lambda)
            Mv = zeros(length(z), size(v, 2));
        end

        function [Mv] = Time_Instance_RHS_zz_Apply(this, v, y, z, t, lambda)
            Mv = zeros(length(z), size(v, 2));
        end

        function [Mv] = Initial_Condition_zz_Apply(this, v, z, lambda)
            A = diag([2 * lambda(1); 6 * lambda(2) * z(2); 12 * lambda(3) * z(3)^2]);
            Mv = A * v;
        end

    end

    methods (Access = public)

        function this = Example_3_Constraint(n_y, n_z, T, n_t)
            this = this@Dynamic_Constraint(n_y, n_z, T, n_t);
        end

    end

end
