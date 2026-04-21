%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Diff_React_Objective < Objective

    properties
        reg_coeff
        m
        M
        T
    end

    methods (Access = public)

        function [val, grad_u, grad_z] = J(this, u, z)
            val = (1 / 2) * (u - this.T)' * this.M * (u - this.T) + (1 / 2) * (this.reg_coeff) * z' * this.M * z;
            grad_u = this.M * (u - this.T);
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

        function this = Diff_React_Objective(diff_react_lofi, reg_coeff)
            this = this@Objective();
            this.reg_coeff = reg_coeff;
            this.m = size(diff_react_lofi.A, 1);
            this.M = diff_react_lofi.pde_meshing.M;
            x = diff_react_lofi.pde_meshing.x;
            y = diff_react_lofi.pde_meshing.y;
            this.T = 15 * exp(-5 * (x.^2 + y.^2)) .* (1 + .5 * sin(pi * x).^2);
        end

    end
end
