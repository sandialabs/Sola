%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Heat_Problem_Objective < Objective
    % Optimal control objective for the heat problem:
    %
    %   min_u  1/2 * \int_0^1 (y - y_d)^2 dx + alpha/2 * \int_0^1 u^2 dx,
    %
    % where y is the state (temperature) and u is the control (distributed
    % heat source).  In Sola notation the state is u and the control is z,
    % so the discretized objective is
    %
    %   J(u, z) = 1/2 (u - y_d)' M (u - y_d) + alpha/2 z' M z.

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

        function this = Heat_Problem_Objective(m, reg_coeff)
            this = this@Objective();
            this.m = m;
            x = linspace(0, 1, m)';
            this.reg_coeff = reg_coeff;

            % Target temperature profile y_d(x), compatible with the
            % homogeneous Dirichlet boundary conditions y(0) = y(1) = 0.
            this.T = sin(pi * x);

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
