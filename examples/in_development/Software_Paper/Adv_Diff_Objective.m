classdef Adv_Diff_Objective < Objective
    % Target-matching objective,

    properties
        m  % State dimension, number of spatial points.
        reg_coeff  % Regularization constant :math:`\beta`.
        T  % Target function T(x) = 50 - 60(x - 1/2)^2
        M  % Mass matrix
    end

    methods (Access = public)

        function [val, grad_u, grad_z] = J(this, u, z)
            grad_u = this.M * (u - this.T);
            grad_z = (this.reg_coeff) * this.M * z;
            val = (1 / 2) * (u - this.T)' * grad_u + (1 / 2) * z' * grad_z;
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

        function this = Adv_Diff_Objective(m, reg_coeff, target)
            arguments
                m
                reg_coeff
                target = []
            end
            this = this@Objective();
            this.m = m;
            x = linspace(0, 1, m)';

            if size(target) == 0
                target = 50 - 60 * (x - 0.5).^2;
            end

            if ~isequal(size(target), size(x))
                error('target and x not aligned');
            end
            this.T = target;
            this.reg_coeff = reg_coeff;
            h = x(2) - x(1);

            M = diag(4 * ones(1, m)) + diag(ones(1, m - 1), 1) + diag(ones(1, m - 1), -1);
            M(1, 1) = .5 * M(1, 1);
            M(end, end) = .5 * M(end, end);
            M = (1 / 6) * h * M;
            this.M = M;
        end

    end
end
