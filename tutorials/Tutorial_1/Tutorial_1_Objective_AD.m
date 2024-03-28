classdef Tutorial_1_Objective_AD < Objective_AD

    properties
        a   % alpha constants.
    end

    methods (Access = public)

        % Constructor: set the alpha constants.
        function this = Tutorial_1_Objective_AD(a, n_u, n_z)
            this = this@Objective_AD(n_u, n_z);
            this.a = a;
        end

        function [val] = J_AD(this, u, z)
            val = sum((u - this.a(1:3)).^2);
            val = val + sum((z - this.a(4:5)).^2);
            val = val + (u(1) * z(1) - this.a(1) * this.a(4))^2;
        end

    end
end
