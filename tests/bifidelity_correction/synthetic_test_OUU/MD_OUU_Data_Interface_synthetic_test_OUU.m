%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef MD_OUU_Data_Interface_synthetic_test_OUU < MD_OUU_Data_Interface

    properties

    end

    %% Pure virtual functions for user implementation
    methods (Access = public)

        function [us_opt] = Load_Optimal_us(this, s)
            us_opt = [];
        end

        function [z_opt] = Load_Optimal_z(this)
            z_opt = load('Optimization_Results.mat', 'z_opt').z_opt;
        end

        function [Xi] = Load_Xi(this)
            Xi = load('Optimization_Results.mat', 'Xi').Xi;
        end

        function [Z] = Load_Z_Data(this)
            Z = [];
        end

        function [Ds] = Load_ds_Data(this, s)
            Ds = [];
        end

    end

    %% Constructor and helper functions
    methods

        function this = MD_OUU_Data_Interface_synthetic_test_OUU()
            this@MD_OUU_Data_Interface(30);
        end

    end

end
