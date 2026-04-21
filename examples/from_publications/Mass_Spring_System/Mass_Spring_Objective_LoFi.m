%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Mass_Spring_Objective_LoFi < Dynamic_Objective

    properties
        reg_coeff
        P_z
    end

    methods

        function this = Mass_Spring_Objective_LoFi(T, n_t)
            n_y = 2;
            n_z = 1;
            this = this@Dynamic_Objective(n_y, n_z, T, n_t);

            this.reg_coeff = 1.e-6;

            P_z = eye(n_t);
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

        function [val] = target(this, t)
            val = 5 * t.^2;
        end

        function [val, grad_y] = g(this, y, t)
            val = 0.5 * (y(1) - this.target(t)).^2;
            grad_y = zeros(2, 1);
            grad_y(1) = y(1) - this.target(t);
        end

        function [val, grad_z] = R(this, z)
            val = 0.5 * this.reg_coeff * (this.w' * (this.P_z * z).^2);
            grad_z = this.reg_coeff * this.P_z' * (this.w .* (this.P_z * z));
        end

        function [Mv] = g_yy_Apply(this, v, y, t)
            Mv = zeros(size(v));
            Mv(1, :) = v(1, :);
        end

        function [Mv] = R_zz_Apply(this, v, z)
            Mv = this.reg_coeff * this.P_z' * diag(this.w) * this.P_z * v;
        end

    end

end
