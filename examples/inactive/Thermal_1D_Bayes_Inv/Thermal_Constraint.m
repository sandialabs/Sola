%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Thermal_Constraint < Constraint

    properties
        m
        x
        dirichlet_bc
        M
        S
        forcing
        B
    end

    methods (Access = public)

        function [u] = State_Solve(this, z)
            D = this.Assembly(z);
            f = this.B * this.forcing + this.dirichlet_bc;
            u = linsolve(D, f);
        end

        function [c] = c(this, u, z)
            D = this.Assembly(z);
            f = this.B * this.forcing + this.dirichlet_bc;
            c = D * u - f;
        end

        function [Mv] = c_u_Transpose_Inverse_Apply(this, v, u, z)
            D = this.Assembly(z);
            Mv = linsolve(D', v);
        end

        function [Mv] = c_z_Transpose_Apply(this, v, u, z)
            D_diff = this.Assembly_z_Jacobian(u);
            Mv = D_diff' * v;
        end

        function [Mv] = c_u_Inverse_Apply(this, v, u, z)
            D = this.Assembly(z);
            Mv = linsolve(D, v);
        end

        function [Mv] = c_z_Apply(this, v, u, z)
            D_diff = this.Assembly_z_Jacobian(u);
            Mv = D_diff * v;
        end

        function [Mv] = c_uu_Apply(this, v, u, z, lambda)
            Mv = zeros(this.m, 1);
        end

        function [Mv] = c_uz_Apply(this, v, u, z, lambda)
            Mv = zeros(this.m, size(v, 2));
            for k = 1:size(v, 2)
                D = this.Assembly(z);
                D_pert = this.Assembly(z + v(:, k));
                Mv(:, k) = D_pert' * lambda - D' * lambda;
            end
        end

        function [Mv] = c_zu_Apply(this, v, u, z, lambda)
            Mv = zeros(this.m, size(v, 2));
            for k = 1:size(v, 2)
                D_diff = this.Assembly_z_Jacobian(u);
                D_diff_pert = this.Assembly_z_Jacobian(u + v(:, k));
                Mv(:, k) = D_diff_pert' * lambda - D_diff' * lambda;
            end
        end

        function [Mv] = c_zz_Apply(this, v, u, z, lambda)
            Mv = zeros(this.m, 1);
        end

    end

    methods (Access = public)

        function [kappa] = Diffusion_Coeff(this, x, z)
            kappa = interp1(this.x, z, x);
        end

        function [D] = Assembly(this, z)
            D = zeros(this.m, this.m);
            h = this.x(2) - this.x(1);
            phi_down_dot = -[1; 1] / h;
            phi_up_dot = [1; 1] / h;
            for i = 1:this.m
                if i > 1
                    x1 = (i - 2) * h + (h / 2) * (-1 / sqrt(3) + 1);
                    x2 = (i - 2) * h + (h / 2) * (1 / sqrt(3) + 1);
                    perm = [this.Diffusion_Coeff(x1, z); this.Diffusion_Coeff(x2, z)];
                    D(i - 1, i) = (h / 2) * sum(phi_up_dot .* phi_down_dot .* perm);
                    D(i, i) = (h / 2) * sum(phi_up_dot .* phi_up_dot .* perm);
                end
                if i < this.m
                    x1 = (i - 1) * h + (h / 2) * (-1 / sqrt(3) + 1);
                    x2 = (i - 1) * h + (h / 2) * (1 / sqrt(3) + 1);
                    perm = [this.Diffusion_Coeff(x1, z); this.Diffusion_Coeff(x2, z)];
                    D(i, i) = D(i, i) + (h / 2) * sum(phi_down_dot .* phi_down_dot .* perm);
                    D(i + 1, i) = (h / 2) * sum(phi_up_dot .* phi_down_dot .* perm);
                end
            end
            D(1, :) = 0 * D(1, :);
            D(this.m, :) = 0 * D(this.m, :);
            D(1, 1) = 1;
            D(this.m, this.m) = 1;
        end

        function [D_diff] = Assembly_z_Jacobian(this, u)
            D_diff = zeros(this.m, this.m);
            h = this.x(2) - this.x(1);
            x1 = (h / 2) * (-1 / sqrt(3) + 1);
            x2 = (h / 2) * (1 / sqrt(3) + 1);
            phi_down = [x2; x1] / h;
            phi_up = [x1; x2] / h;
            phi_down_dot = -[1; 1] / h;
            phi_up_dot = [1; 1] / h;
            for i = 1:this.m
                if i > 1
                    u_prime = u(i - 1) * phi_down_dot + u(i) * phi_up_dot;
                    D_diff(i - 1, i) = (h / 2) * sum(phi_up .* phi_down_dot .* u_prime);
                    D_diff(i, i) = (h / 2) * sum(phi_up .* phi_up_dot .* u_prime);
                end
                if i < this.m
                    u_prime = u(i) * phi_down_dot + u(i + 1) * phi_up_dot;
                    D_diff(i, i) = D_diff(i, i) + (h / 2) * sum(phi_down .* phi_down_dot .* u_prime);
                    D_diff(i + 1, i) = (h / 2) * sum(phi_down .* phi_up_dot .* u_prime);
                end
            end
            D_diff(1, :) = 0 * D_diff(1, :);
            D_diff(this.m, :) = 0 * D_diff(this.m, :);
        end

        function this = Thermal_Constraint(m)
            this = this@Constraint();

            this.m = m;
            this.x = linspace(0, 1, m)';

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

            this.dirichlet_bc = zeros(m, 1);
            this.dirichlet_bc(1) = 0;
            this.dirichlet_bc(m) = 0;

            this.forcing = this.x.^2;

            B = (10^3) * this.M;
            B(1, :) = 0 * B(1, :);
            B(end, :) = 0 * B(end, :);
            this.B = B;

        end

        function [diffs] = Jacobian_Test(this)
            z = randn(this.m, 1);
            u = rand(this.m, 1);
            h = 10.^(-2:-1:-6);
            p = length(h);

            D_diff = this.Assembly_z_Jacobian(u);
            v = randn(this.m, 1);
            D_diff_v = D_diff * v;
            D = this.Assembly(z);

            diffs = zeros(p, 1);
            for k = 1:p
                D_pert = this.Assembly(z + h(k) * v);
                FD = (D_pert * u - D * u) / h(k);
                diffs(k) = norm(D_diff_v - FD);
            end
            disp('Jacobian finite difference check');
            for k = 1:p
                disp(['h = ', num2str(h(k)), ' and error = ', num2str(diffs(k))]);
            end
            disp(' ');
        end

    end
end
