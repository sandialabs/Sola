%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Example_1_Constraint < Parameterized_Constraint

    % Constraint:
    % c(u,z,theta) = [ theta(1) * u(1) + u(2) - theta(2) * z(1)
    %                  z(1) * u(2) - theta(3) * z(2)
    %                  theta(4) * u(3)^3 - z(2)^2               ]

    methods (Access = public)

        function this = Example_1_Constraint(theta)
            this@Parameterized_Constraint(theta);
        end

        function [con] = Parameterized_c(this, u, z, theta)
            con = [theta(1) * u(1) + u(2) - theta(2) * z(1)
                   z(1) * u(2) - theta(3) * z(2)
                   theta(4) * u(3)^3 - z(2)^2];
        end

        function [u] = Parameterized_State_Solve(this, z, theta)
            u2 = theta(3) * z(2) / z(1);
            u1 = (1 / theta(1)) * (theta(2) * z(1) - u2);
            u3 = (z(2)^2 / theta(4))^(1 / 3);
            u = [u1; u2; u3];
        end

        function [u_out] = Parameterized_c_u_Transpose_Inverse_Apply(this, u_in, u, z, theta)
            u_out = [u_in(1, :) / theta(1)
                     (u_in(2, :) - u_in(1, :) / theta(1)) / z(1)
                     u_in(3, :) / (3 * theta(4) * u(3)^2)];
        end

        function [z_out] = Parameterized_c_z_Transpose_Apply(this, u_in, u, z, theta)
            z_out = [-theta(2) * u_in(1, :) + u_in(2, :) * u(2)
                     -theta(3) * u_in(2, :) - 2 * u_in(3, :) * z(2)];
        end

        function [u_out] = Parameterized_c_u_Inverse_Apply(this, u_in, u, z, theta)
            u2z = u_in(2, :) / z(1);
            u_out = [(u_in(1, :) - u2z) / theta(1)
                     u2z
                     u_in(3, :) / (3 * theta(4) * u(3)^2)];
        end

        function [u_out] = Parameterized_c_z_Apply(this, z_in, u, z, theta)
            u_out = [-theta(2) * z_in(1, :)
                     z_in(1, :) * u(2) - z_in(2, :) * theta(3)
                     -2 * z_in(2, :) * z(2)];
        end

        function [u_out] = Parameterized_c_theta_Apply(this, theta_in, u, z, theta)
            u_out = [u(1) * theta_in(1) - z(1) * theta_in(2)
                     -z(2) * theta_in(3)
                     u(3)^3 * theta_in(4)];
        end

        function [u_out] = Parameterized_c_uu_Apply(this, u_in, u, z, lambda, theta)
            u_out = [0; 0; 6 * theta(4) * lambda(3) * u(3) * u_in(3)];
        end

        function [u_out] = Parameterized_c_uz_Apply(this, z_in, u, z, lambda, theta)
            u_out = [0; lambda(2) * z_in(1); 0];
        end

        function [u_out] = Parameterized_c_utheta_Apply(this, theta_in, u, z, lambda, theta)
            u_out = [lambda(1) * theta_in(1)
                     0.0
                     3 * u(3)^2 * lambda(3) * theta_in(4)];
        end

        function [z_out] = Parameterized_c_zu_Apply(this, u_in, u, z, lambda, theta)
            z_out = [lambda(2) * u_in(2); 0];
        end

        function [z_out] = Parameterized_c_ztheta_Apply(this, theta_in, u, z, lambda, theta)
            z_out = [-lambda(1) * theta_in(2)
                     -lambda(2) * theta_in(3)];
        end

        function [z_out] = Parameterized_c_zz_Apply(this, z_in, u, z, lambda, theta)
            z_out = [0; -2 * lambda(3) * z_in(2)];
        end

    end
end
