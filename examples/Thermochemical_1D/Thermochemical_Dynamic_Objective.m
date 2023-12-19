classdef Thermochemical_Dynamic_Objective < Dynamic_Objective

    properties
        con
        reg_coeff
        gl_weights
    end

    methods (Access = public)

        function [val, grad_y] = Time_Instance_Objective(this, y, t)
            grad_y = -this.con.I_u1' * this.con.fe.M * (this.con.I_u1 * y);
            val = (1 / 2) * y' * grad_y;
        end

        function [val, grad_z] = Regularization_Objective(this, z)
            zt = reshape(z, this.con.fe.m, this.con.control_time_nodes);
            val = this.reg_coeff * (1 / 2) * diag(zt' * this.con.fe.M * zt)' * this.gl_weights;
            tmp = (this.con.fe.M * zt) .* (ones(this.con.fe.m, 1) * this.gl_weights');
            grad_z = this.reg_coeff * tmp(:);
        end

        function [Mv] = Time_Instance_Objective_yy_Apply(this, v, y, t)
            Mv = zeros(length(y), size(v, 2));
            for k = 1:size(v, 2)
                [~, Mv(:, k)] = this.Time_Instance_Objective(v(:, k), t);
            end
        end

        function [Mv] = Regularization_Objective_zz_Apply(this, v, z)
            Mv = zeros(length(z), size(v, 2));
            for k = 1:size(v, 2)
                [~, Mv(:, k)] = this.Regularization_Objective(v(:, k));
            end
        end

        function this = Thermochemical_Dynamic_Objective(con, reg_coeff)
            this@Dynamic_Objective(con.n_y, con.n_z, con.T, con.n_t);
            this.con = con;
            this.reg_coeff = reg_coeff;
            this.gl_weights = (1 / 2) * (1 ./ ((0:(this.con.control_time_nodes - 1))' + 1 / 2));
        end

    end

end
