classdef MD_Data_Interface_Py < MD_Data_Interface

    properties
        data_interface_py
    end

    methods

        function [u_opt] = Load_Optimal_u(this)
            u_opt = this.data_interface_py.Load_Optimal_u();
            u_opt = double(u_opt);
            if size(u_opt, 1) == 1
                u_opt = u_opt';
            end
        end

        function [z_opt] = Load_Optimal_z(this)
            z_opt = this.data_interface_py.Load_Optimal_z();
            z_opt = double(z_opt);
            if size(z_opt, 1) == 1
                z_opt = z_opt';
            end
        end

        function [Z] = Load_Z_Data(this)
            Z = this.data_interface_py.Load_Z_Data();
            Z = double(Z);
            if size(Z, 1) == 1
                Z = Z';
            end
        end

        function [D] = Load_d_Data(this)
            D = this.data_interface_py.Load_d_Data();
            D = double(D);
            if size(D, 1) == 1
                D = D';
            end
        end

        function this = MD_Data_Interface_Py(data_interface_py)
            this.data_interface_py = data_interface_py;
        end

    end

end
