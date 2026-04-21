%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Example_1_Constraint < Constraint

    % Algebraic constraint of the form
    % u_1 + u_2     + 0     = z_1
    %
    % 0   + z_1u_2 + 0     = z_2
    %
    % 0   +  0      + u_3^3 = z_2^2

    methods (Access = public)

        function this = Example_1_Constraint()

        end

        function [u] = State_Solve(this, z)
            u = zeros(3, 1);
            u(1) = z(1) - z(2) / z(1);
            u(2) = z(2) / z(1);
            u(3) = z(2).^(2 / 3);
        end

        function [Mv] = c_u_Transpose_Inverse_Apply(this, v, u, z)
            A = [1 1 0; 0 z(1) 0; 0 0 3 * u(3).^2];
            Mv = linsolve(A', v);
        end

        function [Mv] = c_z_Transpose_Apply(this, v, u, z)
            A = [-1, u(2), 0; 0, -1, -2 * z(2)];
            Mv = A * v;
        end

        function [Mv] = c_u_Inverse_Apply(this, v, u, z)
            A = [1, 1, 0; 0, z(1), 0; 0, 0, 3 * u(3)^2];
            Mv = linsolve(A, v);
        end

        function [Mv] = c_z_Apply(this, v, u, z)
            A = [-1, u(2), 0; 0, -1, -2 * z(2)];
            Mv = A' * v;
        end

        function [Mv] = c_uu_Apply(this, v, u, z, lambda)
            A = zeros(3, 3);
            A(3, 3) = 6 * lambda(3) * u(3);
            Mv = A * v;
        end

        function [Mv] = c_uz_Apply(this, v, u, z, lambda)
            A = zeros(length(u), length(z));
            A(2, 1) = lambda(2);
            Mv = A * v;
        end

        function [Mv] = c_zu_Apply(this, v, u, z, lambda)
            A = zeros(length(z), length(u));
            A(1, 2) = lambda(2);
            Mv = A * v;
        end

        function [Mv] = c_zz_Apply(this, v, u, z, lambda)
            A = zeros(2, 2);
            A(2, 2) = -2 * lambda(3);
            Mv = A * v;
        end

    end

end
