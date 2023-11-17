classdef MD_Data_Interface < handle

    properties
        u_opt
        z_opt
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
        end

        % Defaults to return empty array
        % Overload function to load high-fidelity data
        function [D] = Load_d_Data(this)
            D = [];
        end

        function this = MD_Data_Interface()

        end

        function [] = Load_Data(this)
            this.u_opt = this.Load_Optimal_u();
            this.z_opt = this.Load_Optimal_z();
            this.Z = this.Load_Z_Data();
            this.D = this.Load_d_Data();
        end

    end

end
