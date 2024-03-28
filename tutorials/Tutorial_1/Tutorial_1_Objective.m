classdef Tutorial_1_Objective < Objective

    properties
        a   % alpha constants.
    end

    methods (Access = public)

        % Constructor: set the alpha constants.
        function this = Tutorial_1_Objective(a)
            this.a = a;
        end

        function [val, grad_u, grad_z] = J(this, u, z)
            u1z1_minus_a1a4 = (u(1) * z(1)) - (this.a(1) * this.a(4));

            % Compute the value of J(u, z).
            val = sum((u - this.a(1:3)).^2);
            val = val + sum((z - this.a(4:5)).^2);
            val = val + u1z1_minus_a1a4^2;

            % Compute the u gradient of J.
            grad_u = 2 * (u - this.a(1:3));
            grad_u(1) = grad_u(1) + 2 * u1z1_minus_a1a4 * z(1);

            % Compute the z gradient of J.
            grad_z = 2 * (z - this.a(4:5));
            grad_z(1) = grad_z(1) + 2 * u1z1_minus_a1a4 * u(1);
        end

        function [u_out] = J_uu_Apply(this, u_in, u, z)
            u_out = 2 * u_in;
            u_out(1, :) = u_out(1, :) .* (1 + z(1)^2);
        end

        function [u_out] = J_uz_Apply(this, z_in, u, z)
            u_out = zeros(length(u), size(z_in, 2));
            u_out(1, :) = (4 * u(1) * z(1) - 2 * this.a(1) * this.a(4)) .* z_in(1, :);
        end

        function [z_out] = J_zu_Apply(this, u_in, u, z)
            z_out = zeros(length(z), size(u_in, 2));
            z_out(1, :) = (4 * u(1) * z(1) - 2 * this.a(1) * this.a(4)) .* u_in(1, :);
        end

        function [z_out] = J_zz_Apply(this, u_in, u, z)
            z_out = 2 * u_in;
            z_out(1, :) = z_out(1, :) .* (1 + u(1)^2);
        end

    end
end
