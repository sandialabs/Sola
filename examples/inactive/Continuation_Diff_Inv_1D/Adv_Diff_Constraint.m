%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Adv_Diff_Constraint < Parameterized_Constraint

    properties
        m
        x
        M
        S
        B
        robin_bc
        forcing
        alpha
        num_state_solves
        num_adjoint_solves
    end

    methods (Access = public)

        function [u] = Parameterized_State_Solve(this, z, theta)
            D = this.Diff_Assembly(z);
            Vel = this.Velocity_Assembly(theta);
            f = this.B * this.forcing;
            A = D + Vel + this.alpha * this.robin_bc;
            u = linsolve(A, f);
            this.num_state_solves = this.num_state_solves + 1;
        end

        function [u_out] = Parameterized_c_u_Transpose_Inverse_Apply(this, u_in, u, z, theta)
            D = this.Diff_Assembly(z);
            Vel = this.Velocity_Assembly(theta);
            c_u = D + Vel + this.alpha * this.robin_bc;
            u_out = linsolve(c_u', u_in);
            this.num_adjoint_solves = this.num_adjoint_solves + size(u_in, 2);
        end

        function [z_out] = Parameterized_c_z_Transpose_Apply(this, u_in, u, z, theta)
            c_z = this.Diff_Assembly_z_Jacobian(u);
            z_out = c_z' * u_in;
        end

        function [u_out] = Parameterized_c_u_Inverse_Apply(this, u_in, u, z, theta)
            D = this.Diff_Assembly(z);
            Vel = this.Velocity_Assembly(theta);
            c_u = D + Vel + this.alpha * this.robin_bc;
            u_out = linsolve(c_u, u_in);
            this.num_adjoint_solves = this.num_adjoint_solves + size(u_in, 2);
        end

        function [u_out] = Parameterized_c_z_Apply(this, z_in, u, z, theta)
            c_z = this.Diff_Assembly_z_Jacobian(u);
            u_out = c_z * z_in;
        end

        function [u_out] = Parameterized_c_theta_Apply(this, theta_in, u, z, theta)
            c_theta = this.Velocity_Assembly_theta_Jacobian(u);
            u_out = c_theta * theta_in;
        end

        function [con] = Parameterized_c(this, u, z, theta)
            D = this.Diff_Assembly(z);
            Vel = this.Velocity_Assembly(theta);
            f = this.B * this.forcing;
            A = D + Vel + this.alpha * this.robin_bc;
            con = A * u - f;
        end

        function [u_out] = Parameterized_c_uu_Apply(this, u_in, u, z, lambda, theta)
            u_out = 0 * u_in;
        end

        function [u_out] = Parameterized_c_uz_Apply(this, z_in, u, z, lambda, theta)
            c_uz = this.Diff_Assembly_z_Jacobian(lambda);
            u_out = c_uz * z_in;
        end

        function [u_out] = Parameterized_c_utheta_Apply(this, theta_in, u, z, lambda, theta)
            c_utheta = this.Velocity_Assembly_utheta_Hessian(lambda);
            u_out = c_utheta * theta_in;
        end

        function [z_out] = Parameterized_c_zu_Apply(this, u_in, u, z, lambda, theta)
            c_uz = this.Diff_Assembly_z_Jacobian(lambda);
            z_out = c_uz' * u_in;
        end

        function [z_out] = Parameterized_c_ztheta_Apply(this, theta_in, u, z, lambda, theta)
            c_ztheta = zeros(length(z), length(theta));
            z_out = c_ztheta * theta_in;
        end

        function [z_out] = Parameterized_c_zz_Apply(this, z_in, u, z, lambda, theta)
            z_out = 0 * z_in;
        end

    end

    methods (Access = public)

        function [kappa] = Diffusion_Coeff(this, x, z)
            kappa = interp1(this.x, z, x);
        end

        function [v] = Velocity_Coeff(this, x, theta)
            v = interp1(this.x, theta, x);
        end

        function [D] = Diff_Assembly(this, z)
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
        end

        function [D_diff] = Diff_Assembly_z_Jacobian(this, u)
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
        end

        function [V] = Velocity_Assembly(this, theta)
            V = zeros(this.m, this.m);
            h = this.x(2) - this.x(1);
            phi_down_dot = -[1; 1] / h;
            phi_up_dot = [1; 1] / h;
            x1 = (h / 2) * (-1 / sqrt(3) + 1);
            x2 = (h / 2) * (1 / sqrt(3) + 1);
            phi_down = [x2; x1] / h;
            phi_up = [x1; x2] / h;
            for i = 1:this.m
                if i > 1
                    x1 = (i - 2) * h + (h / 2) * (-1 / sqrt(3) + 1);
                    x2 = (i - 2) * h + (h / 2) * (1 / sqrt(3) + 1);
                    vel = [this.Velocity_Coeff(x1, theta); this.Velocity_Coeff(x2, theta)];
                    V(i - 1, i) = (h / 2) * sum(phi_up_dot .* phi_down .* vel);
                    V(i, i) = (h / 2) * sum(phi_up_dot .* phi_up .* vel);
                end
                if i < this.m
                    x1 = (i - 1) * h + (h / 2) * (-1 / sqrt(3) + 1);
                    x2 = (i - 1) * h + (h / 2) * (1 / sqrt(3) + 1);
                    vel = [this.Velocity_Coeff(x1, theta); this.Velocity_Coeff(x2, theta)];
                    V(i, i) = V(i, i) + (h / 2) * sum(phi_down_dot .* phi_down .* vel);
                    V(i + 1, i) = (h / 2) * sum(phi_down_dot .* phi_up .* vel);
                end
            end
        end

        function [V_diff] = Velocity_Assembly_theta_Jacobian(this, u)
            V_diff = zeros(this.m, this.m);
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
                    V_diff(i - 1, i) = (h / 2) * sum(phi_up .* phi_down .* u_prime);
                    V_diff(i, i) = (h / 2) * sum(phi_up .* phi_up .* u_prime);
                end
                if i < this.m
                    u_prime = u(i) * phi_down_dot + u(i + 1) * phi_up_dot;
                    V_diff(i, i) = V_diff(i, i) + (h / 2) * sum(phi_down .* phi_down .* u_prime);
                    V_diff(i + 1, i) = (h / 2) * sum(phi_down .* phi_up .* u_prime);
                end
            end
        end

        function [V_hess] = Velocity_Assembly_utheta_Hessian(this, lambda)
            V_hess = zeros(this.m, this.m);
            h = this.x(2) - this.x(1);
            x1 = (h / 2) * (-1 / sqrt(3) + 1);
            x2 = (h / 2) * (1 / sqrt(3) + 1);
            phi_down = [x2; x1] / h;
            phi_up = [x1; x2] / h;
            phi_down_dot = -[1; 1] / h;
            phi_up_dot = [1; 1] / h;
            for i = 1:this.m
                if i > 1
                    lam = lambda(i - 1) * phi_down + lambda(i) * phi_up;
                    V_hess(i - 1, i) = (h / 2) * sum(phi_up .* phi_down_dot .* lam);
                    V_hess(i, i) = (h / 2) * sum(phi_up .* phi_up_dot .* lam);
                end
                if i < this.m
                    lam = lambda(i) * phi_down + lambda(i + 1) * phi_up;
                    V_hess(i, i) = V_hess(i, i) + (h / 2) * sum(phi_down .* phi_down_dot .* lam);
                    V_hess(i + 1, i) = (h / 2) * sum(phi_down .* phi_up_dot .* lam);
                end
            end
        end

        function this = Adv_Diff_Constraint(theta)
            this = this@Parameterized_Constraint(theta);
            this.m = length(theta);
            m = this.m;
            this.x = linspace(0, 1, m)';
            h = this.x(2) - this.x(1);

            % Mass matrix
            M = diag(4 * ones(1, m)) + diag(ones(1, m - 1), 1) + diag(ones(1, m - 1), -1);
            M(1, 1) = .5 * M(1, 1);
            M(end, end) = .5 * M(end, end);
            M = (1 / 6) * h * M;
            this.M = M;

            % Stiffness matrix
            S = diag(2 * ones(1, m)) + (-1) * diag(ones(1, m - 1), 1) + (-1) * diag(ones(1, m - 1), -1);
            S(1, 1) = .5 * S(1, 1);
            S(end, end) = .5 * S(end, end);
            S = (1 / h) * S;
            this.S = S;

            this.forcing = exp(-50 * (this.x - 0.5).^2);

            B = this.M;
            this.B = B;

            % Robin boundary condition operator
            robin_bc = zeros(m, m);
            robin_bc(1, 1) = 1;
            robin_bc(end, end) = 1;
            this.robin_bc = robin_bc;

            this.alpha = 1;

            this.num_state_solves = 0;
            this.num_adjoint_solves = 0;
        end

        % Method of manufactured solutions test for solver
        function [] = MMS_Check(this)

            this.alpha = 1;

            syms x;
            y = -(4 / 5) * x.^2 + x + 1;
            kappa = 1 + x;
            vel = 1 + x.^2;
            yp = diff(y, x);
            rhs = -diff(kappa * diff(y, x), x) + vel * diff(y, x);
            subs(yp * kappa, x, 0) - this.alpha * subs(y, x, 0);
            -subs(yp * kappa, x, 1) - this.alpha * subs(y, x, 1);

            mms = -(4 / 5) * this.x.^2 + this.x + 1;
            z = 1 + this.x;
            theta = 1 + this.x.^2;
            this.forcing = double(subs(rhs, this.x));

            u = this.Parameterized_State_Solve(z, theta);
            figure;
            plot(this.x, u, this.x, mms, '--');
        end

    end
end
