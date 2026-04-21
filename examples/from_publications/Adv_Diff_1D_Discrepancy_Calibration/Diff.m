%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Diff < Constrained_Optimization

    properties
        m
        reg_coeff
        reg_mat
        T
        diff_coeff
        robin_coeff
        x
        M
        S
        robin_bc
    end

    methods (Access = public)

        %% Pure virtual functions for gradient computation

        function [val, grad_u, grad_z] = Objective(this, u, z)
            val = (1 / 2) * (u - this.T)' * this.M * (u - this.T) + (1 / 2) * (this.reg_coeff) * z' * this.reg_mat * z;
            grad_u = this.M * (u - this.T);
            grad_z = (this.reg_coeff) * this.reg_mat * z;
        end

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
            this = this@Constrained_Optimization();
            this.m = adv_diff.m;
            this.x = adv_diff.x;
            this.T = adv_diff.T;
            this.reg_coeff = adv_diff.reg_coeff;
            this.diff_coeff = adv_diff.diff_coeff;
            this.robin_coeff = adv_diff.robin_coeff;
            this.M = adv_diff.M;
            this.S = adv_diff.S;
            this.robin_bc = adv_diff.robin_bc;
            this.reg_mat = adv_diff.reg_mat;
        end

    end
end
