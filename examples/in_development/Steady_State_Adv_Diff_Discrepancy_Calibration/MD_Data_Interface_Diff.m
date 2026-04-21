%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef MD_Data_Interface_Diff < MD_Data_Interface

    properties
        xi
    end

    methods

        function [u_opt] = Load_Optimal_u(this)
            u_opt = load(['Optimization_Results_xi_', num2str(this.xi), '.mat']).u_lofi;
        end

        function [z_opt] = Load_Optimal_z(this)
            z_opt = load(['Optimization_Results_xi_', num2str(this.xi), '.mat']).z_lofi;
        end

        function [Z] = Load_Z_Data(this)
            Z = load(['Optimization_Results_xi_', num2str(this.xi), '.mat']).Z;
        end

        function [D] = Load_d_Data(this)
            D = load(['Optimization_Results_xi_', num2str(this.xi), '.mat']).D;
        end

        function this = MD_Data_Interface_Diff(xi)
            this.xi = xi;
        end

    end

end
