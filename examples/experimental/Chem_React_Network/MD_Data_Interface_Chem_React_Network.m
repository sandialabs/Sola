classdef MD_Data_Interface_Chem_React_Network < MD_Data_Interface

    properties

    end

    methods

        function [u_opt] = Load_Optimal_u(this)
            u_opt = load('Optimization_Results.mat').u_opt;
        end

        function [z_opt] = Load_Optimal_z(this)
            z_opt = load('Optimization_Results.mat').z_opt;
        end

        function [Z] = Load_Z_Data(this)
            Z = load('Optimization_Results.mat').z_opt;
        end

        function [D] = Load_d_Data(this)
            D = load('Optimization_Results.mat').u_hifi - this.Load_Optimal_u();
        end

        function this = MD_Data_Interface_Chem_React_Network()

        end

    end

end
