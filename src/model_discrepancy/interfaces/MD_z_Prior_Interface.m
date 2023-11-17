classdef MD_z_Prior_Interface < handle

    properties

    end

    methods (Abstract, Access = public)

        %% Pure virtual functions

        [z_out] = Apply_W_z_Inverse(this, z_in)

    end

    methods

        % Factorize W_z^{-1} = F*F^T, function gives z_out = F*z_in
        % This function must be implemented to enable posterior update sampling
        function [z_out] = Apply_W_z_Inverse_Factor(this, z_in)
            z_out = [];
            disp('Apply_W_z_Inverse_Factor must be implemented to use sampling algorithms');
        end

        % Apply W_z matrix
        % This function must be implemented to enable Hessian GEVP
        function [z_out] = Apply_W_z(this, z_in)
            z_out = [];
            disp('Apply_W_z must be implemented to use Hessian GEVP');
        end

        function this = MD_z_Prior_Interface()

        end

    end

end
