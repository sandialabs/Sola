classdef Diff_Objective < Objective

    properties
        diff
        reg_coeff
        m
        M
        T
    end

    methods (Access = public)

        function [val, grad_u, grad_z] = J(this, u, z)
            d =  u - this.T;
            val = (1 / 2) * d' * this.M * d + (1 / 2) * (this.reg_coeff) * z' * this.M * z;
            grad_u = this.M * d;
            grad_z = (this.reg_coeff) * this.M * z;
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
            Mv = this.reg_coeff * this.M * v;
        end

    end

    methods (Access = public)

        function this = Diff_Objective(diff, reg_coeff)
            this = this@Objective();
            this.diff = diff;
            this.reg_coeff = reg_coeff;
            this.m = size(this.diff.A, 1);
            this.M = this.diff.pde_meshing.M;
            x = this.diff.pde_meshing.x;
            y = this.diff.pde_meshing.y;
            this.T = 10 * exp(-50 * ((x + .3).^2 + 3 * (y - .1).^2)) + 15 * exp(-50 * ((x - .4).^2 + 3 * (y - .3).^2));

        end

    end
end
