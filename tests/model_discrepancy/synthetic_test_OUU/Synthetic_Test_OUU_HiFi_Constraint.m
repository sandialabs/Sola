classdef Synthetic_Test_OUU_HiFi_Constraint < Parametric_Constraint

    properties
        alpha
    end

    methods (Access = public)

        function this = Synthetic_Test_OUU_HiFi_Constraint(theta)
            this@Parametric_Constraint(theta);
            this.alpha = 1.2;
        end

        function [con] = Parametric_c(this, u, z, theta)
            con = u - this.alpha * theta(1) * z.^3 - theta(2);
        end

        function [u] = Parametric_State_Solve(this, z, theta)
            u = this.alpha * theta(1) * z.^3 + theta(2);
        end

        function [u_out] = Parametric_c_u_Transpose_Inverse_Apply(this, u_in, u, z, theta)
            u_out = u_in;
        end

        function [z_out] = Parametric_c_z_Transpose_Apply(this, u_in, u, z, theta)
            z_out = -diag(this.alpha * theta(1) * 3 * z.^2) * u_in;
        end

        function [u_out] = Parametric_c_u_Inverse_Apply(this, u_in, u, z, theta)
            u_out = u_in;
        end

        function [u_out] = Parametric_c_z_Apply(this, z_in, u, z, theta)
            u_out = -diag(this.alpha * theta(1) * 3 * z.^2) * z_in;
        end

        function [u_out] = Parametric_c_uu_Apply(this, u_in, u, z, lambda, theta)
            u_out = 0 * u_in;
        end

        function [u_out] = Parametric_c_uz_Apply(this, z_in, u, z, lambda, theta)
            u_out = 0 * z_in;
        end

        function [z_out] = Parametric_c_zu_Apply(this, u_in, u, z, lambda, theta)
            z_out = 0 * u_in;
        end

        function [z_out] = Parametric_c_zz_Apply(this, z_in, u, z, lambda, theta)
            z_out = -diag(this.alpha * theta(1) * 6 * z .* lambda) * z_in;
        end

    end
end