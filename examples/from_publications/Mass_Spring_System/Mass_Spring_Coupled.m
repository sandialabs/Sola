%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Mass_Spring_Coupled < Dynamic_Constraint

    properties
        m1
        m2
        k1
        k2
        k3
        f2
        P_z
    end

    methods

        function this = Mass_Spring_Coupled(T, N)
            n_y = 4;
            n_z = 1;
            this = this@Dynamic_Constraint(n_y, n_z, T, N);

            this.m1 = 1;
            this.m2 = 10;
            this.k1 = 1;
            this.k2 = 1;
            this.k3 = 1;
            this.f2 = 0;

            P_z = eye(N);
            P_z = P_z(:, 2:end);
            P_z(1, 1) = 1;
            this.P_z = P_z;
        end

        function [w] = Temporal_Weights(this, t)
            w = (this.t_mesh - t) / (this.t_mesh(2) - this.t_mesh(1));
            Im = intersect(find(w <= 0), find(abs(w) <= 1));
            Ip = intersect(find(w > 0), find(abs(w) <= 1));
            w(abs(w) > 1) = 0;
            w(Im) = 1 + w(Im);
            w(Ip) = 1 - w(Ip);
        end

        % ODE system with four states y=(x_1,v_1,x_2,v_2)
        % x_1' = v_1
        % v_1' = (1/m_1)*( k_2*x_2 - (k_1+k_2)*x_1 + f_1(z) )
        % x_2' = v_2
        % v_2' = (1/m_2)*( k_2*x_1 - (k_2+k_3)*x_2 + f_2 )
        function [f, f_y, f_z] = f(this, y, z, t)
            x1 = y(1);
            v1 = y(2);
            x2 = y(3);
            v2 = y(4);

            coeffs = this.Temporal_Weights(t);
            f1 = (this.P_z * z)' * coeffs;

            f = zeros(4, 1);
            f(1) = v1;
            f(2) = (1 / this.m1) * (this.k2 * x2 - (this.k1 + this.k2) * x1 + f1);
            f(3) = v2;
            f(4) = (1 / this.m2) * (this.k2 * x1 - (this.k2 + this.k3) * x2 + this.f2);

            f_y = zeros(4, 4);
            f_y(1, 2) = 1;
            f_y(2, 1) = -(1 / this.m1) * (this.k1 + this.k2);
            f_y(2, 3) = (1 / this.m1) * this.k2;
            f_y(3, 4) = 1;
            f_y(4, 1) = (1 / this.m2) * this.k2;
            f_y(4, 3) = -(1 / this.m2) * (this.k2 + this.k3);

            f_z = zeros(4, size(this.P_z, 2));
            f_z(2, :) = (1 / this.m1) * this.P_z' * coeffs;
        end

        function [h, h_z] = h(this, z)
            h = zeros(4, 1);
            h_z = zeros(4, size(this.P_z, 2));
        end

        function [Mv] = f_yy_Apply(this, v, y, z, t, lambda)
            Mv = 0 * v;
        end

        function [Mv] = f_yz_Apply(this, v, y, z, t, lambda)
            Mv = zeros(4, size(v, 2));
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
