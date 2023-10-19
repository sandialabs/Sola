classdef Example_1_Objective < Objective

    %% Problem description
    % \min_{z \in R^2} J(u,z) = || u - (a_1, a_2, a_3)^T ||^2 + || z - (a_4, a_5)^T ||^2 + (u_1 z_1 - a_1 a_4)^2

    properties
        a;
    end

    methods (Access = public)

        function this = Example_1_Objective(a)
            this.a = a;
        end

        function [val, grad_u, grad_z] = J(this, u, z)
            u1z1_minus_a1a4 = (u(1)*z(1)) - (this.a(1)*this.a(4));

            % Calculate the value of J(u, z).
            val = sum((u - this.a(1:3)).^2);
            val = val + sum((z - this.a(4:5)).^2);
            val = val + u1z1_minus_a1a4^2;

            % Calculate the u gradient of J.
            grad_u = 2*(u - this.a(1:3));
            grad_u(1) = grad_u(1) + 2*u1z1_minus_a1a4*z(1);

            % Calculate the z gradient of J.
            grad_z = 2*(z - this.a(4:5));
            grad_z(1) = grad_z(1) + 2*u1z1_minus_a1a4*u(1);
        end

        function [Mv] = J_uu_Apply(this, v, u, z)
            Mv = 2 * v;
            Mv(1, :) = Mv(1, :) .* (1 + z(1)^2);
        end

        function [Mv] = J_uz_Apply(this, v, u, z)
            Mv = zeros(length(u), size(v, 2));
            Mv(1, :) = (4*u(1)*z(1) - 2*this.a(1)*this.a(4)) .* v(1, :);
        end

        function [Mv] = J_zu_Apply(this, v, u, z)
            Mv = zeros(length(z), size(v, 2));
            Mv(1, :) = (4*u(1)*z(1) - 2*this.a(1)*this.a(4)) .* v(1, :);
        end

        function [Mv] = J_zz_Apply(this, v, u, z)
            Mv = 2 * v;
            Mv(1, :) = Mv(1, :) * (1 + u(1)^2);
        end

    end

end
