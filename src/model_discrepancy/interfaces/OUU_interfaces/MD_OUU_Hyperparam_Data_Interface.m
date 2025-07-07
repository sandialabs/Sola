classdef MD_OUU_Hyperparam_Data_Interface < MD_Data_Interface

    properties
        md_ouu_data_interface
    end

    methods (Access = public)

        function [u_opt] = Load_Optimal_u(this)
            tmp = this.md_ouu_data_interface.Load_Optimal_u();
            tmp = this.md_ouu_data_interface.Reshape_State_to_Mat(tmp);
            u_opt = mean(tmp, 2);
        end

        function [z_opt] = Load_Optimal_z(this)
            z_opt = this.md_ouu_data_interface.Load_Optimal_z();
        end

        function [Z] = Load_Z_Data(this)
            Z = this.md_ouu_data_interface.Load_Z_Data();
        end

        function [D] = Load_d_Data(this)
            tmp = this.md_ouu_data_interface.Load_d_Data();
            m = length(this.u_opt);
            N = size(tmp, 2);
            D  = zeros(m, N);
            for k = 1:N
                D(:, k) = mean(this.md_ouu_data_interface.Reshape_State_to_Mat(tmp(:, k)), 2);
            end
        end

        function [I] = Separate_State_Components(this, i)
            I = this.md_ouu_data_interface.Separate_State_Components_Per_Sample(i);
        end

    end

    %% Constructor and helper functions
    methods

        function this = MD_OUU_Hyperparam_Data_Interface(md_ouu_data_interface)
            this.md_ouu_data_interface = md_ouu_data_interface;
        end

    end

end
