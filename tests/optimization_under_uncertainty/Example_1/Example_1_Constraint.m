%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Example_1_Constraint < Parametric_Constraint

    methods (Access = public)

        function this = Example_1_Constraint(theta)
            this@Parametric_Constraint(theta);
        end

        function [con] = Parametric_c(this, u, z, theta)
            con = zeros(3, 1);
            con(1) = u(1) + u(2) - z(1) - theta(1);
            con(2) = z(1) * u(2) - z(2) - theta(2);
            con(3) = u(3)^3 - z(2)^2 - theta(3);
        end

        function [u] = Parametric_State_Solve(this, z, theta)
            u2 = (z(2) + theta(2)) / z(1);
            u = [z(1) - u2 + theta(1); u2; (z(2)^2 + theta(3))^(1 / 3)];
        end

        function [u_out] = Parametric_c_u_Transpose_Inverse_Apply(this, u_in, u, z, theta)
            u_out = [u_in(1, :)
                     (u_in(2, :) - u_in(1, :)) / z(1)
                     u_in(3, :) / (3 * u(3)^2)];
        end

        function [z_out] = Parametric_c_z_Transpose_Apply(this, u_in, u, z, theta)
            z_out = [-u_in(1, :) + u_in(2, :) * u(2)
                     -u_in(2, :) - 2 * u_in(3, :) * z(2)];
        end

        function [u_out] = Parametric_c_u_Inverse_Apply(this, u_in, u, z, theta)
            u2z = u_in(2, :) / z(1);
            u_out = [u_in(1, :) - u2z
                     u2z
                     u_in(3, :) / (3 * u(3)^2)];
        end

        function [u_out] = Parametric_c_z_Apply(this, z_in, u, z, theta)
            u_out = [-z_in(1, :)
                     z_in(1, :) * u(2) - z_in(2, :)
                     -2 * z_in(2, :) * z(2)];
        end

        function [u_out] = Parametric_c_uu_Apply(this, u_in, u, z, lambda, theta)
            u_out = [0; 0; 6 * lambda(3) * u(3) * u_in(3)];
        end

        function [u_out] = Parametric_c_uz_Apply(this, z_in, u, z, lambda, theta)
            u_out = [0; lambda(2) * z_in(1); 0];
        end

        function [z_out] = Parametric_c_zu_Apply(this, u_in, u, z, lambda, theta)
            z_out = [lambda(2) * u_in(2); 0];
        end

        function [z_out] = Parametric_c_zz_Apply(this, z_in, u, z, lambda, theta)
            z_out = [0; -2 * lambda(3) * z_in(2)];
        end

    end
end
