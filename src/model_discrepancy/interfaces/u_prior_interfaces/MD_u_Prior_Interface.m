classdef MD_u_Prior_Interface < handle

    properties

    end

    %% Pure virtual functions for user implementation
    methods (Abstract, Access = public)

        [u_out] = Apply_M_u(this, u_in)

        [u_out] = Apply_W_u_Plus_scalar_M_u_Inverse(this, u_in, scalar)

        [u_out] = Apply_W_u_Inverse(this, u_in)

    end

    %% Virtual functions for user implementation
    methods

        % Compute samples from a mean zero Gaussian with covariance W_u^{-1}
        function [u_out] = Sample_with_Covariance_W_u_Inverse(this, num_samples)
            u_out = [];
            disp('MD_u_Prior_Interface::Sample_with_Covariance_W_u_Inverse must be implemented to use sampling algorithms');
        end

        % Compute samples from a mean zero Gaussian with covariance (W_u+scalar*M_u)^{-1}
        function [u_out] = Sample_with_Covariance_W_u_Plus_scalar_M_u_Inverse(this, num_samples, scalar)
            u_out = [];
            disp('MD_u_Prior_Interface::Sample_with_Covariance_W_u_Plus_scalar_M_u_Inverse must be implemented to use sampling algorithms');
        end

    end

    %% Constructor
    methods

        function this = MD_u_Prior_Interface()

        end

    end

end
