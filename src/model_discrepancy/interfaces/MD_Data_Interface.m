classdef MD_Data_Interface < handle

    properties
        u_opt
        z_opt
        u_init
        z_init
        Z
        D
    end

    methods (Abstract, Access = public)

        %% Pure virtual functions

        [u_opt] = Load_Optimal_u(this)

        [z_opt] = Load_Optimal_z(this)

    end

    methods

        % Defaults to return empty array
        % Overload function to load high-fidelity data
        function [Z] = Load_Z_Data(this)
            Z = [];
            disp('Load_Z_Data must be implemented for any analyses except for optimal experimental design');
        end

        % Defaults to return empty array
        % Overload function to load high-fidelity data
        function [D] = Load_d_Data(this)
            D = [];
            disp('Load_d_Data must be implemented for any analyses except for optimal experimental design');
        end

        function Set_Z_and_D(this, Z, D)
            this.Z = Z;
            this.D = D;
        end

        function Update_z_opt(this, z_opt)
            this.z_opt = z_opt;
        end

        function this = MD_Data_Interface()

        end

        function [] = Load_Data(this)
            % Load Data from Virtual Functions
            if isempty(this.u_init)
                this.u_init = this.Load_Optimal_u();
            end
            if isempty(this.z_init)
                this.z_init = this.Load_Optimal_z();
            end
            if isempty(this.Z)
                this.Z = this.Load_Z_Data();
            end
            if isempty(this.D)
                this.D = this.Load_d_Data();
            end

            % Default u_opt/z_opt to u_init/z_init
            if isempty(this.u_opt)
                this.u_opt = this.u_init;
            end
            if isempty(this.z_opt)
                this.z_opt = this.z_init;
            end
        end

    end

end
