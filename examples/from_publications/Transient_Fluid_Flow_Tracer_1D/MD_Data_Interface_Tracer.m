%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef MD_Data_Interface_Tracer < MD_Data_Interface
    methods

        function [u_opt] = Load_Optimal_u(this)
            u_opt = this.u_lofi;
        end

        function [z_opt] = Load_Optimal_z(this)
            z_opt = this.z_lofi;
        end

        function [Z] = Load_Z_Data(this)
            Z = this.Z;
        end

        function [D] = Load_d_Data(this)
            D = this.D;
        end

        function this = MD_Data_Interface_Tracer(varargin)
            switch nargin
                case 2
                    [this.u_opt, this.z_opt] = deal(varargin{:});
                otherwise
                    error("Please enter the correct number of inputs into the data interface.");
            end

        end

    end

end
