classdef Example_2_Constraint < Dynamic_Constraint
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

        function [f, f_y, f_z] = Time_Instance_RHS(this, y, z, t)
            f = y;
            f_y = eye(2);
            f_z = zeros(2, 2);
        end

        function [h, h_z] = Initial_Condition(this, z)
            h = z;
            h_z = eye(2);
        end

        function [Mv] = Time_Instance_RHS_yy_Apply(this, v, y, z, t, lambda)
            Mv = 0 * v;
        end

        function [Mv] = Time_Instance_RHS_yz_Apply(this, v, y, z, t, lambda)
            Mv = zeros(this.n_y, size(v, 2));
        end

        function [Mv] = Time_Instance_RHS_zy_Apply(this, v, y, z, t, lambda)
            Mv = zeros(length(z), size(v, 2));
        end

        function [Mv] = Time_Instance_RHS_zz_Apply(this, v, y, z, t, lambda)
            Mv = 0 * v;
        end

        function [Mv] = Initial_Condition_zz_Apply(this, v, z, lambda)
            Mv = 0 * z;
        end

    end

    methods (Access = public)

        function this = Example_2_Constraint(n_y, n_z, T, n_t)
            this = this@Dynamic_Constraint(n_y, n_z, T, n_t);
        end

    end

end
