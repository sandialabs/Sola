classdef MD_z_Prior_Interface < handle

    properties

    end

    %% Pure virtual functions for user implementation
    methods (Abstract, Access = public)

        [z_out] = Apply_M_z(this, z_in)

        [z_out] = Apply_W_z_Inverse(this, z_in)

    end

    %% Virtual functions for user implementation
    methods

        % Compute samples from a mean zero Gaussian with covariance W_z^{-1}
        function [z_out] = Sample_with_Covariance_W_z_Inverse(this, num_samples)
            z_out = [];
            disp('MD_z_Prior_Interface::Sample_with_Covariance_W_z_Inverse must be implemented to use sampling algorithms');
        end

        % Apply W_z matrix
        % This function must be implemented to enable Hessian GEVP
        function [z_out] = Apply_W_z(this, z_in)
            z_out = [];
            disp('MD_z_Prior_Interface::Apply_W_z must be implemented to use the Hessian GEVP');
        end

    end

    %% Constructor
    methods

        function this = MD_z_Prior_Interface()

        end

    end

end
