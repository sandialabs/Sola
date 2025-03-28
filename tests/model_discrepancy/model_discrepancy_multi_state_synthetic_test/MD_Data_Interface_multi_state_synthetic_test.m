classdef MD_Data_Interface_multi_state_synthetic_test < MD_Data_Interface

    properties
        m
        x
    end

    methods (Access = public)

        function [u_opt] = Load_Optimal_u(this)
            u_opt_1 = (1 + this.x).^3;
            u_opt_2 = u_opt_1 + 1;
            u_opt = [u_opt_1; u_opt_2];
        end

        function [z_opt] = Load_Optimal_z(this)
            z_opt = 1 + this.x;
        end

        function [Z] = Load_Z_Data(this)
            Z = zeros(this.m, 2);
            Z(:, 1) = 1 + this.x;
            Z(:, 2) = this.x + this.x.^2;
        end

        function [D] = Load_d_Data(this)
            Z = this.Load_Z_Data();
            D = [(1.1 - 1) * (Z.^3); (1.1^2 - 1) * (Z.^3)];
        end

        function [I] = Separate_State_Components(this, i)
            if i == 1
                I = (1:this.m)';
            elseif i == 2
                I = ((this.m + 1):(2 * this.m))';
            else
                disp('Error in Separate_State_Components');
            end
        end

    end

    methods

        function this = MD_Data_Interface_multi_state_synthetic_test(m)
            this.m = m;
            this.x = linspace(0, 1, m)';
        end

    end

end
