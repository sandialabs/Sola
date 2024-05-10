classdef Tutorial_1_Constraint < Constraint

    methods (Access = public)

        function [u] = State_Solve(this, z)
            u2 = z(2) / z(1);
            u = [z(1) - u2; u2; z(2)^(2 / 3)];
        end

        % Jacobian actions.
        function [u_out] = c_u_Transpose_Inverse_Apply(this, u_in, u, z)
            u_out = [u_in(1, :)
                     (u_in(2, :) - u_in(1, :)) / z(1)
                     u_in(3, :) / (3 * u(3)^2)];
        end

        function [z_out] = c_z_Transpose_Apply(this, u_in, u, z)
            z_out = [-u_in(1, :) + u_in(2, :) * u(2)
                     -u_in(2, :) - 2 * u_in(3, :) * z(2)];
        end

        function [u_out] = c_u_Inverse_Apply(this, u_in, u, z)
            u2z = u_in(2, :) / z(1);
            u_out = [u_in(1, :) - u2z
                     u2z
                     u_in(3, :) / (3 * u(3)^2)];
        end

        function [u_out] = c_z_Apply(this, z_in, u, z)
            u_out = [-z_in(1, :)
                     z_in(1, :) * u(2) - z_in(2, :)
                     -2 * z_in(2, :) * z(2)];
        end

        % Hessian actions. These methods are not required when
        % Reduced_Space_Optimization.Gauss_Newton_Hess=true.
        function [u_out] = c_uu_Apply(this, u_in, u, z, lambda)
            u_out = [0; 0; 6 * lambda(3) * u(3) * u_in(3)];
        end

        function [u_out] = c_uz_Apply(this, z_in, u, z, lambda)
            u_out = [0; lambda(2) * z_in(1); 0];
        end

        function [z_out] = c_zu_Apply(this, u_in, u, z, lambda)
            z_out = [lambda(2) * u_in(2); 0];
        end

        function [z_out] = c_zz_Apply(this, z_in, u, z, lambda)
            z_out = [0; -2 * lambda(3) * z_in(2)];
        end

        % This method is required for finite difference checks.
        function [con] = c(this, u, z)
            con = [u(1) + u(2) - z(1)
                   z(1) * u(2) - z(2)
                   u(3)^3 - z(2)^2];
        end

    end
end
