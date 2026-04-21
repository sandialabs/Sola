%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Chem_React_Network_Objective < Dynamic_Objective

    properties
        con
        target
        target_species
        reg_coeff
    end

    methods

        function this = Chem_React_Network_Objective(T, n_t, con)
            n_y = 9;
            n_z = 1;
            this = this@Dynamic_Objective(n_y, n_z, T, n_t);

            this.con = con;
            this.target = [200 * con.state_scale / (con.nA * con.vol)];
            this.target_species = [5];
            this.reg_coeff = 0;
        end

        function [val, grad_y] = g(this, y, t)
            val = 0;
            grad_y = 0 * y;
            if abs(t - this.T) < 1.e-12
                val = (this.con.nA * this.con.vol / this.con.state_scale) * this.n_t * 0.5 * sum((y(this.target_species) - this.target).^2);
                grad_y(this.target_species) = (this.con.nA * this.con.vol / this.con.state_scale) * this.n_t * (y(this.target_species) - this.target);
            end
        end

        function [val, grad_z] = R(this, z)
            val = 0.5 * this.reg_coeff * z^2;
            grad_z = this.reg_coeff * z;
        end

        function [Mv] = g_yy_Apply(this, v, y, t)
            Mv = zeros(size(v));
            if abs(t - this.T) < 1.e-12
                Mv(this.target_species, :) = (this.con.nA * this.con.vol / this.con.state_scale) * this.n_t * v(this.target_species, :);
            end
        end

        function [Mv] = R_zz_Apply(this, v, z)
            Mv = this.reg_coeff * v;
        end

    end

end
