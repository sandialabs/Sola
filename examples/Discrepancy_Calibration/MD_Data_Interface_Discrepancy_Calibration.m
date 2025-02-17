classdef MD_Data_Interface_Discrepancy_Calibration < MD_Data_Interface

    properties
        z_lofi_in
        u_lofi_in
        Z_in
        D_in
    end

    methods

        function [u_opt] = Load_Optimal_u(this)
            u_opt = this.u_lofi_in;
        end

        function [z_opt] = Load_Optimal_z(this)
            z_opt = this.z_lofi_in;
        end

        function [Z] = Load_Z_Data(this)
            Z = this.Z_in;
        end

        function [D] = Load_d_Data(this)
            D = this.D_in;
        end

        function this = MD_Data_Interface_Discrepancy_Calibration(z_lofi,u_lofi,Z,D)
            this.z_lofi_in = z_lofi;
            this.u_lofi_in = u_lofi;
            this.Z_in = Z;
            this.D_in = D;
        end

    end

end
