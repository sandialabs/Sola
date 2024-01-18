classdef Thermochemical_Data_Interface < MD_Data_Interface

    properties

    end

    methods (Access = public)

        function [u_opt] = Load_Optimal_u(this)
            u_opt = load('LoFi_Opt_Results.mat', 'u').u;
        end

        function [z_opt] = Load_Optimal_z(this)
            z_opt = load('LoFi_Opt_Results.mat', 'z').z;
        end

        function [Z] = Load_Z_Data(this)
            Z = load('Discrepancy_Evaluations.mat', 'Z').Z;
        end

        function [D] = Load_d_Data(this)
            D = load('Discrepancy_Evaluations.mat', 'D').D;
        end

        function this = Thermochemical_Data_Interface()

        end

    end

end
