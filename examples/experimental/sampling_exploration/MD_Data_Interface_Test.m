%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef MD_Data_Interface_Test < MD_Data_Interface

    properties
        m
    end

    %% Pure virtual functions for user implementation
    methods

        function [u_opt] = Load_Optimal_u(this)
            u_opt = ones(this.m,1);
        end

        function [z_opt] = Load_Optimal_z(this)
            z_opt = ones(this.m,1);
        end


        function [Z] = Load_Z_Data(this)
            Z = ones(this.m,1);
        end

        function [D] = Load_d_Data(this)
            D = 0.3 * ones(this.m,1);
        end

    end

    %% Constructor and helper functions
    methods

        function this = MD_Data_Interface_Test(m)
            this.m = m;
        end
    end

end
