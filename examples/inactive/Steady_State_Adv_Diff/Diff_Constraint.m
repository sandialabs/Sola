classdef Diff_Constraint < Constraint

    properties
        m
        diff_coeff
        robin_coeff
        x
        M
        S
        robin_bc
    end

    methods (Access = public)

        function [u] = State_Solve(this, z)
            A = this.diff_coeff * this.S + this.robin_coeff * this.robin_bc;
            b = (10^2) * this.M * z;
            u = linsolve(A, b);
        end

        function [Mv] = c_u_Transpose_Inverse_Apply(this, v, u, z)
            A = this.diff_coeff * this.S + this.robin_coeff * this.robin_bc;
            Mv = linsolve(A', v);
        end

        function [Mv] = c_z_Transpose_Apply(this, v, u, z)
            Mv = -(10^2) * this.M' * v;
        end

        function [Mv] = c_u_Inverse_Apply(this, v, u, z)
            A = this.diff_coeff * this.S + this.robin_coeff * this.robin_bc;
            Mv = linsolve(A, v);
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
            this.x = adv_diff.x;
            this.diff_coeff = adv_diff.diff_coeff;
            this.robin_coeff = adv_diff.robin_coeff;
            this.M = adv_diff.M;
            this.S = adv_diff.S;
            this.robin_bc = adv_diff.robin_bc;
        end

    end
end
