%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Adv_Diff < Constrained_Optimization

    properties
        m
        reg_coeff
        reg_mat
        T
        diff_coeff
        vel_coeff
        robin_coeff
        x
        M
        S
        V
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
            A = this.diff_coeff * this.S + this.vel_coeff * this.V + this.robin_coeff * this.robin_bc;
            b = (10^2) * this.M * z;
            u = linsolve(A, b);
        end

        function [Mv] = c_u_Transpose_Inverse_Apply(this, v, u, z)
            A = this.diff_coeff * this.S + this.vel_coeff * this.V + this.robin_coeff * this.robin_bc;
            Mv = linsolve(A', v);
        end

        function [Mv] = c_z_Transpose_Apply(this, v, u, z)
            Mv = -(10^2) * this.M' * v;
        end

        function [Mv] = c_u_Inverse_Apply(this, v, u, z)
            A = this.diff_coeff * this.S + this.vel_coeff * this.V + this.robin_coeff * this.robin_bc;
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

        function this = Adv_Diff(m, diff_coeff, vel_coeff, robin_coeff, reg_coeff)
            this = this@Constrained_Optimization();
            this.m = m;
            this.diff_coeff = diff_coeff;
            this.vel_coeff = vel_coeff;
            this.robin_coeff = robin_coeff;
            this.x = linspace(0, 1, m)';
            this.reg_coeff = reg_coeff;
            this.T = 50 - 60 * (this.x - 0.5).^2;

            h = this.x(2) - this.x(1);

            M = diag(4 * ones(1, m)) + diag(ones(1, m - 1), 1) + diag(ones(1, m - 1), -1);
            M(1, 1) = .5 * M(1, 1);
            M(end, end) = .5 * M(end, end);
            M = (1 / 6) * h * M;
            this.M = M;

            S = diag(2 * ones(1, m)) + (-1) * diag(ones(1, m - 1), 1) + (-1) * diag(ones(1, m - 1), -1);
            S(1, 1) = .5 * S(1, 1);
            S(end, end) = .5 * S(end, end);
            S = (1 / h) * S;
            this.S = S;

            V = diag(0 * ones(1, m)) + (1 / 2) * diag(ones(1, m - 1), 1) + (-1 / 2) * diag(ones(1, m - 1), -1);
            V(1, 1) = -1 / 2;
            V(end, end) = 1 / 2;
            this.V = V;

            robin_bc = zeros(m, m);
            robin_bc(1, 1) = 1;
            robin_bc(end, end) = 1;
            this.robin_bc = robin_bc;

            this.reg_mat = M + (10^-3) * S;
        end

    end
end
