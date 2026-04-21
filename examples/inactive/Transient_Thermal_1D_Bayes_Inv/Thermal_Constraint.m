%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Thermal_Constraint < Dynamic_Constraint

    properties
        x
        M
        S
        forcing
    end

    methods (Access = public)

        function [f, f_y, f_z] = f(this, y, z, t)
            D = this.Assembly(z);
            f_y = -linsolve(this.M, D);
            f = f_y * y + (10^3) * this.forcing(this.x, t);
            f_z = -linsolve(this.M, this.Assembly_z_Jacobian(y));
        end

        function [h, h_z] = h(this, z)
            h = ones(this.n_y, 1);
            h_z = zeros(this.n_y, this.n_y);
        end

        function [Mv] = f_yy_Apply(this, v, y, z, t, lambda)
            Mv = 0 * v;
        end

        function [Mv] = f_yz_Apply(this, v, y, z, t, lambda)
            Mv = zeros(this.n_y, size(v, 2));
            for k = 1:size(v, 2)
                D = this.Assembly(z);
                D_pert = this.Assembly(z + v(:, k));
                Mv(:, k) = -(D_pert' - D') * linsolve(this.M, lambda);
            end
        end

        function [Mv] = f_zy_Apply(this, v, y, z, t, lambda)
            Mv = zeros(this.n_y, size(v, 2));
            for k = 1:size(v, 2)
                D_diff = this.Assembly_z_Jacobian(y);
                D_diff_pert = this.Assembly_z_Jacobian(y + v(:, k));
                Mv(:, k) = -(D_diff_pert' - D_diff') * linsolve(this.M, lambda);
            end
        end

        function [Mv] = f_zz_Apply(this, v, y, z, t, lambda)
            Mv = 0 * v;
        end

        function [Mv] = h_zz_Apply(this, v, z, lambda)
            Mv = 0 * v;
        end

    end

    methods (Access = public)

        function [val] = Forcing_Function(this, x, t)
            val = this.forcing(x, t);
        end

        function [kappa] = Diffusion_Coeff(this, x, z)
            kappa = interp1(this.x, z, x);
        end

        function [D] = Assembly(this, z)
            h = this.x(2) - this.x(1);
            x1 = (0:(this.n_y - 2)) * h + (h / 2) * (-1 / sqrt(3) + 1);
            x2 = (0:(this.n_y - 2)) * h + (h / 2) * (1 / sqrt(3) + 1);
            diff_x = this.Diffusion_Coeff([x1; x2], z);

            s = sum(diff_x, 1);
            D = diag(([0, s] + [s, 0]) * (1 / h) / 2) + (-1) * diag(s, 1) * (1 / h) / 2 + (-1) * diag(s, -1) * (1 / h) / 2;

            %             phi_down_dot = -[1;1]/h;
            %             phi_up_dot = [1;1]/h;
            %             diff_x1 = diff_x(1:(this.n_y-1));
            %             diff_x2 = diff_x(this.n_y:end);
            %             D = zeros(this.n_y,this.n_y);
            %             for i = 1:this.n_y
            %                 if i > 1
            %                     perm = [diff_x1(i-1) ; diff_x2(i-1)];
            %                     D(i-1,i) = (h/2)*sum(phi_up_dot.*phi_down_dot.*perm);
            %                     D(i,i) = (h/2)*sum(phi_up_dot.*phi_up_dot.*perm);
            %                 end
            %                 if i < this.n_y
            %                     perm = [diff_x1(i) ; diff_x2(i)];
            %                     D(i,i) = D(i,i) + (h/2)*sum(phi_down_dot.*phi_down_dot.*perm);
            %                     D(i+1,i) = (h/2)*sum(phi_up_dot.*phi_down_dot.*perm);
            %                 end
            %             end
        end

        function [D_diff] = Assembly_z_Jacobian(this, u)
            h = this.x(2) - this.x(1);
            up = (u(2:end) - u(1:end - 1)) / h;
            d = [-(1 / 2) * up(1); (1 / 2) * up(1:end - 1) - (1 / 2) * up(2:end); (1 / 2) * up(end)];
            D_diff = diag(d) + (-1 / 2) * diag(up, 1) + (1 / 2) * diag(up, -1);

            %             D_diff = zeros(this.n_y,this.n_y);
            %             x1 = (h/2)*(-1/sqrt(3) + 1);
            %             x2 = (h/2)*(1/sqrt(3) + 1);
            %             phi_down = [x2;x1]/h;
            %             phi_up = [x1;x2]/h;
            %             phi_down_dot = -[1;1]/h;
            %             phi_up_dot = [1;1]/h;
            %             for i = 1:this.n_y
            %                 if i > 1
            %                     u_prime = u(i-1)*phi_down_dot + u(i)*phi_up_dot;
            %                     D_diff(i-1,i) = (h/2)*sum(phi_up.*phi_down_dot.*u_prime);
            %                     D_diff(i,i) = (h/2)*sum(phi_up.*phi_up_dot.*u_prime);
            %                 end
            %                 if i < this.n_y
            %                     u_prime = u(i)*phi_down_dot + u(i+1)*phi_up_dot;
            %                     D_diff(i,i) = D_diff(i,i) + (h/2)*sum(phi_down.*phi_down_dot.*u_prime);
            %                     D_diff(i+1,i) = (h/2)*sum(phi_down.*phi_up_dot.*u_prime);
            %                 end
            %             end
        end

        function this = Thermal_Constraint(n_y, n_z, T, n_t)
            this = this@Dynamic_Constraint(n_y, n_z, T, n_t);

            this.x = linspace(0, 1, n_y)';

            h = this.x(2) - this.x(1);

            M = diag(4 * ones(1, n_y)) + diag(ones(1, n_y - 1), 1) + diag(ones(1, n_y - 1), -1);
            M(1, 1) = .5 * M(1, 1);
            M(end, end) = .5 * M(end, end);
            M = (1 / 6) * h * M;
            this.M = M;

            S = diag(2 * ones(1, n_y)) + (-1) * diag(ones(1, n_y - 1), 1) + (-1) * diag(ones(1, n_y - 1), -1);
            S(1, 1) = .5 * S(1, 1);
            S(end, end) = .5 * S(end, end);
            S = (1 / h) * S;
            this.S = S;

            this.forcing = @(x, t) this.x.^2;

        end

    end
end
