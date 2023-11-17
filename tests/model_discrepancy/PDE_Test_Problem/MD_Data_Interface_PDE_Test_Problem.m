classdef MD_Data_Interface_PDE_Test_Problem < MD_Data_Interface

    properties

    end

    methods (Access = public)

        function [u_opt] = Load_Optimal_u(this)
            u_opt = load('u_opt.mat').u_opt;
        end

        function [z_opt] = Load_Optimal_z(this)
            z_opt = load('z_opt.mat').z_opt;
        end

        function [Z] = Load_Z_Data(this)
            Z = load('Z.mat').Z;
        end

        function [D] = Load_d_Data(this)
            D = load('D.mat').D;
        end

    end

    methods

        function this = MD_Data_Interface_PDE_Test_Problem()

        end

    end

end
