classdef Thermochemical_Data_Interface < MD_Data_Interface

    properties
        n_y
        n_t
    end

    methods (Access = public)

        function [u_opt] = Load_Optimal_u(this)
            u_opt = load('LoFi_Opt_Results.mat', 'u').u;
        end

        function [z_opt] = Load_Optimal_z(this)
            z_opt = load('LoFi_Opt_Results.mat', 'z').z;
        end

        function [Z] = Load_Z_Data(this)
            Z = load('HiFi_Opt_Results.mat', 'z').z - load('LoFi_Opt_Results.mat', 'z').z;
        end

        function [D] = Load_d_Data(this)
            u_hifi = load('HiFi_Opt_Results.mat', 'u').u;
            u_lofi = load('LoFi_Opt_Results.mat', 'u').u;
            u_hifi_rs = reshape(u_hifi, 6 * this.n_y, this.n_t);
            u_hifi_rs = u_hifi_rs(kron([1; 2; 4; 5], (1:this.n_y)'), :);
            D = u_hifi_rs(:) - u_lofi;
        end

        function this = Thermochemical_Data_Interface(n_y, n_t)
            this.n_y = n_y;
            this.n_t = n_t;
        end

    end

end
