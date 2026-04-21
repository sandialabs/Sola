%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef Thermal_Objective < Dynamic_Objective

    methods (Access = public)

        function [val, grad_y] = g(this, y, t)
            val = (1 / 2) * (y' * y);
            grad_y = y;
        end

        function [val, grad_z] = R(this, z)
            val = (1 / 2) * (z' * z);
            grad_z = z;
        end

        function [Mv] = g_yy_Apply(this, v, y, t)
            Mv = v;
        end

        function [Mv] = R_zz_Apply(this, v, z)
            Mv = v;
        end

    end

    methods (Access = public)

        function this = Thermal_Objective(m, n, T, N)
            this@Dynamic_Objective(m, n, T, N);
        end

    end

end
