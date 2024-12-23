classdef Chem_React_Network_Objective < Dynamic_Objective

    properties
        con
        target
        reg_coeff
    end

    methods

        function this = Chem_React_Network_Objective(T, n_t,con)
            n_y = 9;
            n_z = 1;
            this = this@Dynamic_Objective(n_y, n_z, T, n_t);

            this.con = con;
            this.target = 125*con.state_scale/(con.nA*con.vol);
            this.reg_coeff = 0;
        end

        function [val, grad_y] = g(this, y, t)
            val = 0;
            grad_y = 0 * y;
            if abs(t-this.T) < 1.e-12
                val = (this.con.nA*this.con.vol/this.con.state_scale) * this.n_t * 0.5 * (y(6)-this.target)^2;
                grad_y(6) = (this.con.nA*this.con.vol/this.con.state_scale) * this.n_t * (y(6)-this.target);
            end
        end

        function [val, grad_z] = R(this, z)
            val = 0.5 * this.reg_coeff * z^2;
            grad_z = this.reg_coeff * z;
        end

        function [Mv] = g_yy_Apply(this, v, y, t)
            Mv = zeros(size(v));
            if abs(t-this.T) < 1.e-12
                Mv(6,:) = (this.con.nA*this.con.vol/this.con.state_scale) * this.n_t * v(6,:);
            end
        end

        function [Mv] = R_zz_Apply(this, v, z)
            Mv = this.reg_coeff * v;
        end

    end

end
