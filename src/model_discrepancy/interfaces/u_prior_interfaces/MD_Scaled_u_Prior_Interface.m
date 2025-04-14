classdef MD_Scaled_u_Prior_Interface < MD_u_Prior_Interface

    %%%%%%%%%%%%%%%%% Express covariance matrix as W_u^{-1} = alpha_u * \acute{W}_u^{-1} so that alpha_u may be isolated %%%%%%%%%%%%%%%%%

    properties
        alpha_u
    end

    %% Pure virtual functions for user implementation
    methods (Abstract, Access = public)

        [u_out] = Apply_M_u(this, u_in)

        [u_out] = Apply_W_u_Acute_Plus_scalar_M_u_Inverse(this, u_in, scalar)

        [u_out] = Apply_W_u_Acute_Inverse(this, u_in)

    end

    %% Virtual functions for user implementation
    methods

        % Compute samples from a mean zero Gaussian with covariance \acute{W}_u^{-1}
        function [u_out] = Sample_with_Covariance_W_u_Acute_Inverse(this, num_samples)
            u_out = [];
            disp('MD_u_Prior_Interface::Sample_with_Covariance_W_u_Acute_Inverse must be implemented to use sampling algorithms');
        end

        % Compute samples from a mean zero Gaussian with covariance (\acute{W}_u+scalar*M_u)^{-1}
        function [u_out] = Sample_with_Covariance_W_u_Acute_Plus_scalar_M_u_Inverse(this, num_samples, scalar)
            u_out = [];
            disp('MD_u_Prior_Interface::Sample_with_Covariance_W_u_Acute_Plus_scalar_M_u_Inverse must be implemented to use sampling algorithms');
        end

    end

    %% Constructor and helper functions
    methods

        function this = MD_Scaled_u_Prior_Interface(alpha_u)
            arguments
                alpha_u (1,1) {mustBeNumeric}
            end
            this.alpha_u = alpha_u;
        end

        function [] = Set_alpha_u(this, alpha_u_new)
            this.alpha_u = alpha_u_new;
        end

    end

    %% Implementation of base class functions
    methods

        function [u_out] = Apply_W_u_Inverse(this, u_in)
            u_out = this.alpha_u * this.Apply_W_u_Acute_Inverse(u_in);
        end

        function [u_out] = Apply_W_u_Plus_scalar_M_u_Inverse(this, u_in, scalar)
            u_out = this.alpha_u * this.Apply_W_u_Acute_Plus_scalar_M_u_Inverse(u_in, scalar * this.alpha_u);
        end

        function [u_out] = Sample_with_Covariance_W_u_Inverse(this, num_samples)
            u_out = sqrt(this.alpha_u) * this.Sample_with_Covariance_W_u_Acute_Inverse(num_samples);
        end

        function [u_out] = Sample_with_Covariance_W_u_Plus_scalar_M_u_Inverse(this, num_samples, scalar)
            u_out = sqrt(this.alpha_u) * this.Sample_with_Covariance_W_u_Acute_Plus_scalar_M_u_Inverse(num_samples, this.alpha_u * scalar);
        end

    end

end
