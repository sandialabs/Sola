%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Diff < Parametric_Constraint

    properties
        m
        diff_coeff
        x
        M
        S
        A
        B
        c
    end

    methods (Access = public)

        %% Pure virtual functions for gradient computation

        function [u] = Parametric_State_Solve(this, z, theta)
            if norm(theta - this.theta_current) ~= 0
                this.Assemble(theta);
            end
            u = linsolve(this.A, this.B * z + this.c);
        end

        function [Mv] = Parametric_c_u_Transpose_Inverse_Apply(this, v, u, z, theta)
            if norm(theta - this.theta_current) ~= 0
                this.Assemble(theta);
            end
            Mv = linsolve(this.A', v);
        end

        function [Mv] = Parametric_c_z_Transpose_Apply(this, v, u, z, theta)
            Mv = -this.B' * v;
        end

        function [Mv] = Parametric_c_u_Inverse_Apply(this, v, u, z, theta)
            if norm(theta - this.theta_current) ~= 0
                this.Assemble(theta);
            end
            Mv = linsolve(this.A, v);
        end

        function [Mv] = Parametric_c_z_Apply(this, v, u, z, theta)
            Mv = -this.B * v;
        end

        function [Mv] = Parametric_c_uu_Apply(this, v, u, z, lambda, theta)
            Mv = zeros(this.m, 1);
        end

        function [Mv] = Parametric_c_uz_Apply(this, v, u, z, lambda, theta)
            Mv = zeros(this.m, 1);
        end

        function [Mv] = Parametric_c_zu_Apply(this, v, u, z, lambda, theta)
            Mv = zeros(this.m, 1);
        end

        function [Mv] = Parametric_c_zz_Apply(this, v, u, z, lambda, theta)
            Mv = zeros(this.m, 1);
        end

    end

    methods (Access = public)

        function [] = Assemble(this, theta)
            this.diff_coeff = theta;
            A = this.diff_coeff * this.S;
            A(1, :) = 0;
            A(1, 1) = 1;
            this.A = A;
        end

        function this = Diff(adv_diff)
            this = this@Parametric_Constraint(adv_diff.theta_current);
            this.m = adv_diff.m;
            this.x = adv_diff.x;
            this.M = adv_diff.M;
            this.S = adv_diff.S;

            B = (10^2) * this.M;
            B(1, :) = 0;
            this.B = B;

            c = zeros(this.m, 1);
            c(1) = 1;
            this.c = c;

            this.Assemble(adv_diff.theta_current);
        end

    end
end
