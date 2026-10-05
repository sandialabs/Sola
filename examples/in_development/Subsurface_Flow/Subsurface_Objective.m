%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Subsurface_Objective < Objective
    % Pressure-targeting objective for the subsurface flow example:
    %
    %   min_q  1/2 * int_0^1 (p - p_d)^2 dx
    %        + beta/2 * int_0^1 q^2 dx.
    %
    % In Sola notation the state is u (pressure p) and the control is z
    % (source/sink rate q), so the discretized objective is
    %
    %   J(u,z) = 1/2 (u - p_d)' M (u - p_d) + beta/2 z' M z.

    properties
        m
        reg_coeff
        M
        reg_mat
        T
    end

    methods (Access = public)

        function [val, grad_u, grad_z] = J(this, u, z)
            val = (1 / 2) * (u - this.T)' * this.M * (u - this.T) ...
                + (1 / 2) * this.reg_coeff * z' * this.reg_mat * z;
            grad_u = this.M * (u - this.T);
            grad_z = this.reg_coeff * this.reg_mat * z;
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

        function this = Subsurface_Objective(m, reg_coeff, p0)
            if nargin < 3
                p0 = 0;
            end

            this = this@Objective();
            this.m = m;
            this.reg_coeff = reg_coeff;

            x = linspace(0, 1, m)';

            % Desired pressure perturbation.  It is compatible with the
            % Dirichlet pressure p = p0 at x = 0 and x = 1.
            this.T = p0 + 0.5 * sin(pi * x) + 0.15 * sin(2 * pi * x);

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
