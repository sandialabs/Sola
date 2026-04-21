%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Diff < Constraint

    properties
        m
        T
        diff_coeff
        xi
        x
        M
        S
        A
        B
        c
    end

    methods (Access = public)

        %% Pure virtual functions for gradient computation

        function [u] = State_Solve(this, z)
            u = linsolve(this.A, this.B * z + this.c);
        end

        function [Mv] = c_u_Transpose_Inverse_Apply(this, v, u, z)
            Mv = linsolve(this.A', v);
        end

        function [Mv] = c_z_Transpose_Apply(this, v, u, z)
            Mv = -this.B' * v;
        end

        function [Mv] = c_u_Inverse_Apply(this, v, u, z)
            Mv = linsolve(this.A, v);
        end

        function [Mv] = c_z_Apply(this, v, u, z)
            Mv = -this.B * v;
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

        function [Mv] = Objective_uu_Apply(this, v, u, z)
            Mv = this.M * v;
        end

        function [Mv] = Objective_uz_Apply(this, v, u, z)
            Mv = zeros(this.m, 1);
        end

        function [Mv] = Objective_zu_Apply(this, v, u, z)
            Mv = zeros(this.m, 1);
        end

        function [Mv] = Objective_zz_Apply(this, v, u, z)
            Mv = this.reg_coeff * this.reg_mat * v;
        end

    end

    methods (Access = public)

        function this = Diff(adv_diff)
            this = this@Constraint();
            this.m = adv_diff.m;
            this.xi = adv_diff.xi;
            this.x = adv_diff.x;
            this.T = adv_diff.T;
            this.diff_coeff = adv_diff.diff_coeff;
            this.M = adv_diff.M;
            this.S = adv_diff.S;

            A = this.diff_coeff * this.S;
            A(1, :) = 0;
            A(1, 1) = 1;
            this.A = A;

            B = (10^2) * this.M;
            B(1, :) = 0;
            this.B = B;

            c = zeros(this.m, 1);
            c(1) = 1;
            this.c = c;
        end

    end
end
