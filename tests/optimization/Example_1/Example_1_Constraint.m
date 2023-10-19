classdef Example_1_Constraint < Constraint

    methods (Access = public)

        function [u] = State_Solve(this, z)
            u2 = z(2) / z(1);
            u = [ z(1) - u2; u2; z(2)^(2/3) ];
        end

        function [Mv] = c_u_Transpose_Inverse_Apply(this, v, u, z)
            Mv = [ v(1, :);
                   (v(2, :) - v(1, :))/z(1);
                   v(3, :) / (3 * u(3)^2) ];
        end

        function [Mv] = c_z_Transpose_Apply(this, lambda, u, z)
            Mv = [ -lambda(1, :) + lambda(2, :) * u(2);
                   -lambda(2, :) - 2*lambda(3, :) * z(2) ];
        end

        function [Mv] = c_u_Inverse_Apply(this, lambda, u, z)
            Mv2 = lambda(2, :) / z(1);
            Mv = [ lambda(1, :) - Mv2;
                   Mv2;
                   lambda(3, :) / (3 * u(3)^2) ];
        end

        function [Mv] = c_z_Apply(this, v, u, z)
            Mv = [ -v(1, :);
                   v(1, :) * u(2) - v(2, :);
                   -2*v(2, :) * z(2) ];
        end

        function [Mv] = c_uu_Apply(this, lambda, u, z, v)
            Mv = [ 0; 0; 6 * lambda(3) * u(3) * v(3)];
        end

        function [Mv] = c_uz_Apply(this, lambda, u, z, v)
            Mv = [ 0; lambda(2)*v(1); 0];
        end

        function [Mv] = c_zu_Apply(this, lambda, u, z, v)
            Mv = [ 0; lambda(1) * v(2) ];
        end

        function [Mv] = c_zz_Apply(this, lambda, u, z, v)
            Mv = [ 0; -2 * lambda(2) * v(2) ];
        end
    end
end
