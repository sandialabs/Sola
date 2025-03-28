classdef MD_Data_Interface < handle

    properties
        u_opt
        z_opt
        Z
        D
        data_shift
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

        % Defaults to return all state elements
        % Overload function to extract component i from the state
        % Returns a vector of integers I index elements of component i
        function [I] = Separate_State_Components(this, i)
            I = 1:length(this.u_opt);
        end

        function this = MD_Data_Interface()

        end

        function [] = Load_Data(this)
            this.u_opt = this.Load_Optimal_u();
            this.z_opt = this.Load_Optimal_z();
            this.Z = this.Load_Z_Data();
            this.D = this.Load_d_Data();
            this.data_shift = zeros(size(this.D, 1), 1);
        end

    end

end
