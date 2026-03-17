classdef Diff_React_Objective < Objective

    properties
        m
        reg_coeff
        M
        reg_mat
        T
    end

    methods (Access = public)

        %% Pure virtual functions for gradient computation

        function [val, grad_u, grad_z] = J(this, u, z)
            val = (1 / 2) * (u - this.T)' * this.M * (u - this.T) + (1 / 2) * (this.reg_coeff) * z' * this.reg_mat * z;
            grad_u = this.M * (u - this.T);
            grad_z = (this.reg_coeff) * this.reg_mat * z;
        end

        function [Mv] = J_uu_Apply(this, v, u, z)
            Mv = this.M * v;
        end

        function [Mv] = J_uz_Apply(this, v, u, z)
            Mv = zeros(this.m, 1);
        end

        function [Mv] = J_zu_Apply(this, v, u, z)
            Mv = zeros(this.m, 1);
        end

        function [Mv] = J_zz_Apply(this, v, u, z)
            Mv = this.reg_coeff * this.reg_mat * v;
        end

    end

    methods (Access = public)

        function this = Diff_React_Objective(m, reg_coeff)
            this = this@Objective();
            this.m = m;
            x = linspace(0, 1, m)';
            this.reg_coeff = reg_coeff;
            this.T = 20 * (x + .5) .* (1.3 - x);

            h = x(2) - x(1);

            M = diag(4 * ones(1, m)) + diag(ones(1, m - 1), 1) + diag(ones(1, m - 1), -1);
            M(1, 1) = .5 * M(1, 1);
            M(end, end) = .5 * M(end, end);
            M = (1 / 6) * h * M;
            this.M = M;

            this.reg_mat = M;
        end

    end
end
