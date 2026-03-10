classdef MD_Data_Interface_Mass_Spring < MD_Data_Interface

    properties

    end

    methods

        function [u_opt] = Load_Optimal_u(this)
            u_opt = load('Optimization_Results.mat').u_lofi;
        end

        function [z_opt] = Load_Optimal_z(this)
            z_opt = load('Optimization_Results.mat').z_tilde;
        end

        function [Z] = Load_Z_Data(this)
            Z = load('Optimization_Results.mat').Z;
        end

        function [D] = Load_d_Data(this)
            D = load('Optimization_Results.mat').D;
        end

        function this = MD_Data_Interface_Mass_Spring()

        end

    end

end
