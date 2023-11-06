classdef Thermal_Objective_AD < Dynamic_Objective_AD

    properties

    end

    methods (Access = public)

        function [val] = Time_Instance_Objective_AD(this, y, t)
            val = (1 / 2) * (y' * y);
        end

        function [val] = Regularization_Objective_AD(this, z)
            val = (1 / 2) * (z' * z);
        end

    end

    methods (Access = public)

        function this = Thermal_Objective_AD(m, n, T, N)
            this@Dynamic_Objective_AD(m, n, T, N);
        end

    end

end
