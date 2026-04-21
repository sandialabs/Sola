%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef MD_OUU_Hyperparam_Data_Interface < MD_Data_Interface

    properties
        ouu_data_interface
    end

    methods (Access = public)

        function [u_opt] = Load_Optimal_u(this)
            tmp = this.ouu_data_interface.Load_Optimal_u();
            tmp = this.ouu_data_interface.Reshape_State_to_Mat(tmp);
            u_opt = mean(tmp, 2);
        end

        function [z_opt] = Load_Optimal_z(this)
            z_opt = this.ouu_data_interface.Load_Optimal_z();
        end

        function [Z] = Load_Z_Data(this)
            Z = this.ouu_data_interface.Load_Z_Data();
        end

        function [D] = Load_d_Data(this)
            tmp = this.ouu_data_interface.Load_d_Data();
            m = length(this.u_opt);
            N = size(tmp, 2);
            D  = zeros(m, N);
            for k = 1:N
                D(:, k) = mean(this.ouu_data_interface.Reshape_State_to_Mat(tmp(:, k)), 2);
            end
        end

        function [I] = Separate_State_Components(this, i)
            I = this.ouu_data_interface.Separate_State_Components_Per_Sample(i);
        end

    end

    %% Constructor and helper functions
    methods

        function this = MD_OUU_Hyperparam_Data_Interface(ouu_data_interface)
            arguments
                ouu_data_interface MD_OUU_Data_Interface
            end
            this.ouu_data_interface = ouu_data_interface;
        end

    end

end
