classdef Thermal_Objective < Dynamic_Objective

    methods (Access = public)

        function [val, grad_y] = Time_Instance_Objective(this, y, t)
            val = (1 / 2) * (y' * y);
            grad_y = y;
        end

        function [val, grad_z] = Regularization_Objective(this, z)
            val = (1 / 2) * (z' * z);
            grad_z = z;
        end

        function [Mv] = Time_Instance_Objective_yy_Apply(this, v, y, t)
            Mv = v;
        end

        function [Mv] = Regularization_Objective_zz_Apply(this, v, z)
            Mv = v;
        end

    end

    methods (Access = public)

        function this = Thermal_Objective(m, n, T, N)
            this@Dynamic_Objective(m, n, T, N);
        end

    end

end
