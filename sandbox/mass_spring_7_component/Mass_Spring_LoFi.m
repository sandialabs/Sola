classdef Mass_Spring_LoFi < Dynamic_Constraint

    properties
        mass_spring_coupled
        m1
        m2
        F
        k
    end

    methods

        function this = Mass_Spring_LoFi(mass_spring_coupled)
            this = this@Dynamic_Constraint(4, mass_spring_coupled.n_z, mass_spring_coupled.T, mass_spring_coupled.n_t);
            this.mass_spring_coupled = mass_spring_coupled;

            this.m1 = mass_spring_coupled.m(1);
            this.m2 = sum(mass_spring_coupled.m(2:end));
            this.F = sum(mass_spring_coupled.F);
            this.k = mass_spring_coupled.k(1);
        end

        function [f, f_y, f_z] = f(this, y, z, t)
            x1 = y(1);
            v1 = y(2);
            x2 = y(3);
            v2 = y(4);

            coeffs = this.mass_spring_coupled.Temporal_Weights(t);
            Fz = (this.mass_spring_coupled.P_z * z)' * coeffs;

            f = zeros(4, 1);
            f(1) = v1;
            f(2) = (1 / this.m1) * ( this.k*(x2-x1) - Fz );
            f(3) = v2;
            f(4) = (1 / this.m2) * ( this.k*(x1-x2) + this.F );

            f_y = zeros(4, 4);
            f_y(1, 2) = 1;
            f_y(2, 1) = (1 / this.m1) * (-this.k);
            f_y(2, 3) = (1 / this.m1) * (this.k);
            f_y(3,4) = 1;
            f_y(4,1) = (1 / this.m2) * (this.k);
            f_y(4,3) = (1 / this.m2) * (-this.k);

            f_z = zeros(4, length(z));
            f_z(2, :) = -(1 / this.m1) * this.mass_spring_coupled.P_z' * coeffs;
        end

        function [h, h_z] = h(this, z)
            h = zeros(4, 1);
            h_z = zeros(4, length(z));
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
