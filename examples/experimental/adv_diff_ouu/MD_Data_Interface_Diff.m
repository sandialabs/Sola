classdef MD_Data_Interface_Diff < MD_OUU_Data_Interface

    properties

    end

    methods

        function [us_opt] = Load_Optimal_us(this, s)
            us_opt = load('Optimization_Results.mat', 'u_lofi').u_lofi(:, s);
        end

        function [z_opt] = Load_Optimal_z(this)
            z_opt = load('Optimization_Results.mat', 'z_lofi').z_lofi;
        end

        function [Z] = Load_Z_Data(this)
            Z = load('Optimization_Results.mat', 'Z').Z;
        end

        function [Ds] = Load_ds_Data(this, s)
            Ds = load('Optimization_Results.mat', 'D').D(:, s, :);
            N = size(Ds, 3);
            Ds = reshape(Ds, [], N);
        end

        function this = MD_Data_Interface_Diff()
            this@MD_OUU_Data_Interface(3);
        end

    end

end
