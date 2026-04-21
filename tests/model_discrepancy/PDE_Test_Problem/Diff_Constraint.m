%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Diff_Constraint < Constraint

    properties
        m
        diff_coeff
        robin_coeff
        x
        M
        S
        robin_bc
        z_lofi
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

        function this = Diff_Constraint(adv_diff_obj, adv_diff_con)
            this = this@Constraint();
            this.m = adv_diff_con.m;
            this.x = adv_diff_con.x;
            this.diff_coeff = adv_diff_con.diff_coeff;
            this.robin_coeff = adv_diff_con.robin_coeff;
            this.M = adv_diff_con.M;
            this.S = adv_diff_con.S;
            this.robin_bc = adv_diff_con.robin_bc;

            A = linsolve(this.diff_coeff * this.S + this.robin_coeff * this.robin_bc, (10^2) * this.M);
            z_lofi = linsolve(A' * this.M * A + adv_diff_obj.reg_coeff * adv_diff_obj.reg_mat, A' * this.M * adv_diff_obj.T);
            this.z_lofi = z_lofi;
        end

    end
end
