classdef Adv_Diff_Objective < Dynamic_Objective

    properties
        M
        x
        beta_reg
    end

    methods (Access = public)

        function [val, grad_y] = g(this, y, t)
            target = this.Evaluate_Target(t, this.x);
            val = .5 * (y - target)' * this.M * (y - target);
            grad_y = this.M * (y - target);
        end

        function [val, grad_z] = R(this, z)
            grad_z = this.beta_reg * this.M * z;
            val = .5 * z' * grad_z;
        end

        function [Mv] = g_yy_Apply(this, v, y, t)
            Mv = this.M * v;
        end

        function [Mv] = R_zz_Apply(this, v, z)
            Mv = this.beta_reg * this.M * v;
        end

    end

    methods (Access = public)

        function this = Adv_Diff_Objective(m, n, T, N)
            this = this@Dynamic_Objective(m, n, T, N);

            % Spatial domain
            this.x = linspace(0, 1, m)';
            h = this.x(2) - this.x(1);

            M = diag(4 * ones(1, m)) + diag(ones(1, m - 1), 1) + diag(ones(1, m - 1), -1);
            M(1, 1) = .5 * M(1, 1);
            M(end, end) = .5 * M(end, end);
            M = (1 / 6) * h * M;

            this.M = M;
            this.beta_reg = 10^-3;

        end

        function [target] = Evaluate_Target(this, t, x)
            target = 0.2 * t^2 * exp(-10 * (x - .5).^2);
        end

    end

end
