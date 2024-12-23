classdef React_Rate_Eqn < Dynamic_Constraint

    properties
        k1;
        k2;
        k3;
        k4;
        k5;
        k6;
        vol;
        nA;
        state_scale;
        opt_scale;
    end

    methods

        function this = React_Rate_Eqn(T,n_t)
            n_y = 9;
            n_z = 1;
            this = this@Dynamic_Constraint(n_y, n_z, T, n_t);

            this.k1 = 1.e9;
            this.k2 = 1.e9;
            this.k3 = 1.e8;
            this.k4 = 2.e8;
            this.k5 = 1.e9;
            this.k6 = 1.e8;
            this.vol = 1.e-11;
            this.nA = 6.023e23;
            this.state_scale = 1.e10;
            this.opt_scale = 1.e-3;
        end

        function [f, f_y, f_z] = f(this, y, z, t)

            y = y/this.state_scale;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            f = zeros(9,1);

            f(1) = -this.k1*y(1)*y(2) - this.k3*y(1)*y(3);
            f(2) = -this.k1*y(1)*y(2) - this.k4*y(2)*y(4);
            f(3) = -this.k2*y(3)*y(4) - this.k3*y(1)*y(3);
            f(4) = -this.k2*y(3)*y(4) - this.k4*y(2)*y(4);

            f(5) = this.k3*y(1)*y(3) - this.k5*y(5)*y(6);
            f(6) = this.k4*y(2)*y(4) - this.k5*y(5)*y(6);

            f(7) =  this.k1*y(1)*y(2) - this.k6*y(7)*y(8);
            f(8) =  this.k2*y(3)*y(4) - this.k6*y(7)*y(8);

            f(9) =  this.k5*y(5)*y(6) + this.k6*y(7)*y(8);

            f = f*this.state_scale;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            f_y = zeros(9,9);

            f_y(1,1) = -this.k1*y(2) - this.k3*y(3);
            f_y(1,2) = -this.k1*y(1);
            f_y(1,3) = -this.k3*y(1);
            f_y(2,1) = -this.k1*y(2);
            f_y(2,2) = -this.k1*y(1) - this.k4*y(4);
            f_y(2,4) = -this.k4*y(2);
            f_y(3,1) = -this.k3*y(3);
            f_y(3,3) = -this.k2*y(4) - this.k3*y(1);
            f_y(3,4) = -this.k2*y(3);
            f_y(4,2) = -this.k4*y(4);
            f_y(4,3) = -this.k2*y(4);
            f_y(4,4) = -this.k2*y(3) - this.k4*y(2);

            f_y(5,1) = this.k3*y(3);
            f_y(5,3) = this.k3*y(1);
            f_y(5,5) = -this.k5*y(6);
            f_y(5,6) = -this.k5*y(5);
            f_y(6,2) = this.k4*y(4);
            f_y(6,4) = this.k4*y(2);
            f_y(6,5) = -this.k5*y(6);
            f_y(6,6) = -this.k5*y(5);

            f_y(7,1) = this.k1*y(2);
            f_y(7,2) = this.k1*y(1);
            f_y(7,7) = -this.k6*y(8);
            f_y(7,8) = -this.k6*y(7);
            f_y(8,3) = this.k2*y(4);
            f_y(8,4) = this.k2*y(3);
            f_y(8,7) = -this.k6*y(8);
            f_y(8,8) = -this.k6*y(7);

            f_y(9,5) = this.k5*y(6);
            f_y(9,6) = this.k5*y(5);
            f_y(9,7) = this.k6*y(8);
            f_y(9,8) = this.k6*y(7);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            f_z = zeros(9,1);
        end

        function [h, h_z] = h(this, z)
            h = zeros(9,1);
            h(1) = (z/this.opt_scale)/(this.nA*this.vol);
            h(2) = 1200/(this.nA*this.vol);
            h(3) = 1200/(this.nA*this.vol);
            h(4) = 1200/(this.nA*this.vol);

            h = h*this.state_scale;

            h_z = zeros(9,1);
            h_z(1) = 1/(this.nA*this.vol);

            h_z = h_z*this.state_scale/this.opt_scale;
        end

        function [Mv] = f_yy_Apply(this, v, y, z, t, lambda)
            M = zeros(9,9);
            M(1,2) = -lambda(1) * this.k1 - lambda(2) * this.k1 + lambda(7) * this.k1;
            M(1,3) = -lambda(1) * this.k3 - lambda(3) * this.k3 + lambda(5) * this.k3;
            M(2,4) = -lambda(2) * this.k4 - lambda(4) * this.k4 + lambda(6) * this.k4;
            M(3,4) = -lambda(3) * this.k2 - lambda(4) * this.k2 + lambda(8) * this.k2;
            M(5,6) = -lambda(5) * this.k5 - lambda(6) * this.k5 + lambda(9) * this.k5;
            M(7,8) = -lambda(7) * this.k6 - lambda(8) * this.k6 + lambda(9) * this.k6;

            M = M + M';
            Mv = M * v;

            Mv = Mv/this.state_scale;
        end

        function [Mv] = f_yz_Apply(this, v, y, z, t, lambda)
            Mv = zeros(9,1);
        end

        function [Mv] = f_zy_Apply(this, v, y, z, t, lambda)
            Mv = 0;
        end

        function [Mv] = f_zz_Apply(this, v, y, z, t, lambda)
            Mv = 0;
        end

        function [Mv] = h_zz_Apply(this, v, z, lambda)
            Mv = 0;
        end

    end

end
