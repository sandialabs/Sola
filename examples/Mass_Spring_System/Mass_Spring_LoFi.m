classdef Mass_Spring_LoFi < Dynamic_Constraint

    properties
        mass_spring_coupled
        m1
        m2
        k1
        k2
        k3
        f2
        P_z
    end

    methods

        function this = Mass_Spring_LoFi(mass_spring_coupled)
            this = this@Dynamic_Constraint(2, mass_spring_coupled.n_z, mass_spring_coupled.T, mass_spring_coupled.n_t);
            this.mass_spring_coupled = mass_spring_coupled;

            this.m1 = mass_spring_coupled.m1;
            this.m2 = mass_spring_coupled.m2;
            this.k1 = mass_spring_coupled.k1;
            this.k2 = mass_spring_coupled.k2;
            this.k3 = mass_spring_coupled.k3;
            this.f2 = mass_spring_coupled.f2;
            this.P_z = mass_spring_coupled.P_z;
        end

        % ODE system with four states y=(x_1,v_1)
        % x_1' = v_1
        % v_1' = (1/m_1)*( k_2*x_2 - (k_1+k_2)*x_1 + f_1(z) )
        function [f, f_y, f_z] = Time_Instance_RHS(this, y, z, t)
            x1 = y(1);
            v1 = y(2);

            coeffs = this.mass_spring_coupled.Temporal_Weights(t);
            f1 = (this.P_z * z)' * coeffs;

            f = zeros(2, 1);
            f(1) = v1;
            f(2) = (1 / this.m1) * (-(this.k1 + this.k2) * x1 + f1);

            f_y = zeros(2, 2);
            f_y(1, 2) = 1;
            f_y(2, 1) = -(1 / this.m1) * (this.k1 + this.k2);

            f_z = zeros(2, size(this.P_z, 2));
            f_z(2, :) = (1 / this.m1) * this.P_z' * coeffs;
        end

        function [h, h_z] = Initial_Condition(this, z)
            h = zeros(2, 1);
            h_z = zeros(2, size(this.P_z, 2));
        end

        function [Mv] = Time_Instance_RHS_yy_Apply(this, v, y, z, t, lambda)
            Mv = 0 * v;
        end

        function [Mv] = Time_Instance_RHS_yz_Apply(this, v, y, z, t, lambda)
            Mv = zeros(2, size(v, 2));
        end

        function [Mv] = Time_Instance_RHS_zy_Apply(this, v, y, z, t, lambda)
            Mv = zeros(length(z), size(v, 2));
        end

        function [Mv] = Time_Instance_RHS_zz_Apply(this, v, y, z, t, lambda)
            Mv = 0 * v;
        end

        function [Mv] = Initial_Condition_zz_Apply(this, v, z, lambda)
            Mv = 0 * v;
        end

    end

end
