classdef Diff_Constraint < Constraint

    properties
        m
        M
        S
        A
    end

    methods (Access = public)

        function [u] = State_Solve(this, z)
            u = linsolve(this.A, (10^2) * this.M * z);
        end

        function [Mv] = c_u_Transpose_Inverse_Apply(this, v, u, z)
            Mv = linsolve(this.A', v);
        end

        function [Mv] = c_z_Transpose_Apply(this, v, u, z)
            Mv = -(10^2) * this.M' * v;
        end

        function [Mv] = c_u_Inverse_Apply(this, v, u, z)
            Mv = linsolve(this.A, v);
        end

        function [Mv] = c_z_Apply(this, v, u, z)
            Mv = -(10^2) * this.M * v;
        end

        function [Mv] = c_uu_Apply(this, v, u, z, lambda)
            Mv = zeros(this.m, 1);
        end

        function [Mv] = c_uz_Apply(this, v, u, z, lambda)
            Mv = zeros(this.m, 1);
        end

        function [Mv] = c_zu_Apply(this, v, u, z, lambda)
            Mv = zeros(this.m, 1);
        end

        function [Mv] = c_zz_Apply(this, v, u, z, lambda)
            Mv = zeros(this.m, 1);
        end

    end

    methods (Access = public)

        function this = Diff_Constraint(adv_diff)
            this = this@Constraint();
            this.m = adv_diff.m;
            this.M = adv_diff.M;
            this.S = adv_diff.S;
            this.A = adv_diff.diff_coeff * adv_diff.S + adv_diff.robin_coeff * adv_diff.robin_bc;
        end

    end
end
