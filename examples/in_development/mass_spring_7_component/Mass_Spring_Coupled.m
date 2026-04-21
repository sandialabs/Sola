%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Mass_Spring_Coupled < Dynamic_Constraint

    properties
        P_z
        m
        k
        F
    end

    methods

        function this = Mass_Spring_Coupled(T, N)
            n_y = 14;
            n_z = 1;
            this = this@Dynamic_Constraint(n_y, n_z, T, N);

            P_z = eye(N);
            P_z = P_z(:, 2:end);
            P_z(1, 1) = 1;
            this.P_z = P_z;

            this.m = zeros(7, 1);
            this.m(1) = 1;
            this.m(2) = 1;
            this.m(3) = 1;
            this.m(4) = 1;
            this.m(5) = .5;
            this.m(6) = .5;
            this.m(7) = .5;

            this.k = zeros(6, 1);
            this.k(1) = 1;
            this.k(2) = 10;
            this.k(3) = 10;
            this.k(4) = 1;
            this.k(5) = 1;
            this.k(6) = 1;

            this.F = zeros(3, 1);
            this.F(1) = 1;
            this.F(2) = 1;
            this.F(3) = 1;
        end

        function [w] = Temporal_Weights(this, t)
            w = (this.t_mesh - t) / (this.t_mesh(2) - this.t_mesh(1));
            Im = intersect(find(w <= 0), find(abs(w) <= 1));
            Ip = intersect(find(w > 0), find(abs(w) <= 1));
            w(abs(w) > 1) = 0;
            w(Im) = 1 + w(Im);
            w(Ip) = 1 - w(Ip);
        end

        % ODE system with seven masses -> 14 states y=(x_1,v_1,x_2,v_2,...)
        % x_i' = v_i
        % v_i' = (1/m_i)*F_i
        function [f, f_y, f_z] = f(this, y, z, t)

            coeffs = this.Temporal_Weights(t);
            Fz = (this.P_z * z)' * coeffs;

            x = y(1:2:end);
            v = y(2:2:end);

            fx = v;
            fv = zeros(7, 1);
            fv(1) = (1 / this.m(1)) * (this.k(1) * (x(2) - x(1)) - Fz);
            fv(2) = (1 / this.m(2)) * (this.k(1) * (x(1) - x(2)) + this.k(2) * (x(3) - x(2)) + this.k(3) * (x(4) - x(2)) + this.k(4) * (x(5) - x(2)));
            fv(3) = (1 / this.m(3)) * (this.k(2) * (x(2) - x(3)) + this.k(5) * (x(6) - x(3)));
            fv(4) = (1 / this.m(4)) * (this.k(3) * (x(2) - x(4)) + this.k(6) * (x(7) - x(4)));
            fv(5) = (1 / this.m(5)) * (this.k(4) * (x(2) - x(5)) + this.F(1) * t);
            fv(6) = (1 / this.m(6)) * (this.k(5) * (x(3) - x(6)) + this.F(2) * t);
            fv(7) = (1 / this.m(7)) * (this.k(6) * (x(4) - x(7)) + this.F(3) * t);

            f = [fx, fv]';
            f = f(:);

            f_y = zeros(14, 14);
            f_y(1:2:end, 2:2:end) = eye(7);

            fv_x = zeros(7, 7);

            fv_x(1, 1) = (1 / this.m(1)) * (-this.k(1));
            fv_x(1, 2) = (1 / this.m(1)) * (this.k(1));

            fv_x(2, 1) = (1 / this.m(2)) * this.k(1);
            fv_x(2, 2) = (1 / this.m(2)) * (-this.k(1) - this.k(2) - this.k(3) - this.k(4));
            fv_x(2, 3) = (1 / this.m(2)) * this.k(2);
            fv_x(2, 4) = (1 / this.m(2)) * this.k(3);
            fv_x(2, 5) = (1 / this.m(2)) * this.k(4);

            fv_x(3, 2) = (1 / this.m(3)) * this.k(2);
            fv_x(3, 3) = (1 / this.m(3)) * (-this.k(2) - this.k(5));
            fv_x(3, 6) = (1 / this.m(3)) * this.k(5);

            fv_x(4, 2) = (1 / this.m(4)) * this.k(3);
            fv_x(4, 4) = (1 / this.m(4)) * (-this.k(3) - this.k(6));
            fv_x(4, 7) = (1 / this.m(4)) * this.k(6);

            fv_x(5, 2) = (1 / this.m(5)) * this.k(4);
            fv_x(5, 5) = (1 / this.m(5)) * (-this.k(4));

            fv_x(6, 3) = (1 / this.m(6)) * this.k(5);
            fv_x(6, 6) = (1 / this.m(6)) * (-this.k(5));

            fv_x(7, 4) = (1 / this.m(7)) * this.k(6);
            fv_x(7, 7) = (1 / this.m(7)) * (-this.k(6));

            f_y(2:2:end, 1:2:end) = fv_x;

            f_z = zeros(14, length(z));
            f_z(2, :) = -(1 / this.m(1)) * this.P_z' * coeffs;

        end

        function [h, h_z] = h(this, z)
            h = zeros(14, 1);
            h_z = zeros(14, length(z));
        end

        function [Mv] = f_yy_Apply(this, v, y, z, t, lambda)
            Mv = 0 * v;
        end

        function [Mv] = f_yz_Apply(this, v, y, z, t, lambda)
            Mv = zeros(14, size(v, 2));
        end

        function [Mv] = f_zy_Apply(this, v, y, z, t, lambda)
            Mv = zeros(length(z), size(v, 2));
        end

        function [Mv] = f_zz_Apply(this, v, y, z, t, lambda)
            Mv = 0 * v;
        end

        function [Mv] = h_zz_Apply(this, v, z, lambda)
            Mv = 0 * v;
        end

    end

end
