%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Synthetic_Test_OUU_Objective < Objective

    properties
        m
        x
        M
        T
    end

    methods (Access = public)

        function this = Synthetic_Test_OUU_Objective(m)
            this.m = m;
            this.x = linspace(0, 1, m)';

            h = this.x(2) - this.x(1);
            M = diag(4 * ones(1, m)) + diag(ones(1, m - 1), 1) + diag(ones(1, m - 1), -1);
            M(1, 1) = .5 * M(1, 1);
            M(end, end) = .5 * M(end, end);
            M = (1 / 6) * h * M;
            this.M = M;

            this.T = (1 + this.x).^3;
        end

        function [val, grad_u, grad_z] = J(this, u, z)
            grad_u = this.M * (u - this.T);
            grad_z = 0 * z;
            val = (1 / 2) * (u - this.T)' * grad_u;
        end

        function [u_out] = J_uu_Apply(this, u_in, u, z)
            u_out = this.M * u_in;
        end

        function [u_out] = J_uz_Apply(this, z_in, u, z)
            u_out = zeros(length(u), size(z_in, 2));
        end

        function [z_out] = J_zu_Apply(this, u_in, u, z)
            z_out = zeros(length(z), size(u_in, 2));
        end

        function [z_out] = J_zz_Apply(this, z_in, u, z)
            z_out = zeros(length(z), size(z_in, 2));
        end

    end

end
