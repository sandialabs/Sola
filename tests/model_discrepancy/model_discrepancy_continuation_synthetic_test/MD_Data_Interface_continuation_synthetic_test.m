classdef MD_Data_Interface_continuation_synthetic_test < MD_Data_Interface

    properties

    end

    methods (Access = public)

        function [u_opt] = Load_Optimal_u(this)
            u_opt = load('Opt_Data.mat', 'u').u;
        end

        function [z_opt] = Load_Optimal_z(this)
            z_opt = load('Opt_Data.mat', 'z').z;
        end

        function [Z] = Load_Z_Data(this)
            Z = load('Opt_Data.mat', 'Z').Z;
        end

        function [D] = Load_d_Data(this)
            D = load('Opt_Data.mat', 'D').D;
        end

    end

    methods

        function this = MD_Data_Interface_continuation_synthetic_test()

        end

    end

end
