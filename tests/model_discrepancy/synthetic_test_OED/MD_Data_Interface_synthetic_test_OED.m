%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef MD_Data_Interface_synthetic_test_OED < MD_Data_Interface

    properties
        m
        x
    end

    methods (Access = public)

        function [u_opt] = Load_Optimal_u(this)
            u_opt = (1 + this.x).^3;
        end

        function [z_opt] = Load_Optimal_z(this)
            z_opt = 1 + this.x;
        end

        function [D] = Evaluate_Discrepancy(this, Z)
            D = .2 * (Z.^3);
        end

        function [Z] = Load_Z_Data(this)
            % Define dummy function to avoid warning
            Z = [];
        end

        function [D] = Load_d_Data(this)
            % Define dummy function to avoid warning
            D = [];
        end

    end

    methods

        function this = MD_Data_Interface_synthetic_test_OED(m)
            this.m = m;
            this.x = linspace(0, 1, m)';
        end

    end

end
