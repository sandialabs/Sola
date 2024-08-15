classdef MD_Data_Interface_Tracer < MD_Data_Interface
    methods

        function [u_opt] = Load_Optimal_u(this)
            u_opt = load('Optimization_Results.mat').u_lofi;
        end

        function [z_opt] = Load_Optimal_z(this)
            z_opt = load('Optimization_Results.mat').z_lofi;
        end

        function [Z] = Load_Z_Data(this)
            Z = load('Optimization_Results.mat').Z;
        end

        function [D] = Load_d_Data(this)
            D = load('Optimization_Results.mat').D;
        end

        function this = MD_Data_Interface_Tracer(varargin)
            switch nargin
                case 2
                    [this.u_init, this.z_init] = deal(varargin{:});
                otherwise
                    error("Please enter the correct number of inputs into the data interface.");
            end

            if isempty(this.z_opt)
                this.z_opt = this.z_init;
            end
        end

    end

end
