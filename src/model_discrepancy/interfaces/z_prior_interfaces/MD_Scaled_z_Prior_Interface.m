classdef MD_Scaled_z_Prior_Interface < MD_z_Prior_Interface

    %%%%%%%%%%%%%%%%% Express covariance matrix as W_z^{-1} = alpha_z * \acute{W}_z^{-1} so that alpha_z may be isolated %%%%%%%%%%%%%%%%%

    properties
        alpha_z
    end

    %% Pure virtual functions for user implementation
    methods (Abstract, Access = public)

        [z_out] = Apply_M_z(this, z_in)

        [z_out] = Apply_W_z_Acute_Inverse(this, z_in)

    end

    %% Virtual functions for user implementation
    methods (Access = public)

        function [z_out] = Sample_with_Covariance_W_z_Acute_Inverse(this, num_samples)
            z_out = [];
            disp('MD_z_Prior_Interface::Sample_with_Covariance_W_z_Acute_Inverse must be implemented to use sampling algorithms');
        end

        function [z_out] = Apply_W_z_Acute(this, z_in)
            z_out = [];
            disp('MD_z_Prior_Interface::Apply_W_z_Acute must be implemented to use the Hessian GEVP');
        end

    end

    %% Constructor and helper functions
    methods

        function this = MD_Scaled_z_Prior_Interface(alpha_z)
            arguments
                alpha_z (1,1) {mustBeNumeric}
            end
            this.alpha_z = alpha_z;
        end

        function [] = Set_alpha_z(this, alpha_z_new)
            this.alpha_z = alpha_z_new;
        end

    end

    %% Implementation of base class virtual functions
    methods

        function [z_out] = Apply_W_z_Inverse(this, z_in)
            z_out = this.alpha_z * this.Apply_W_z_Acute_Inverse(z_in);
        end

        function [z_out] = Sample_with_Covariance_W_z_Inverse(this, num_samples)
            z_out = this.Sample_with_Covariance_W_z_Acute_Inverse(num_samples);
            z_out = sqrt(this.alpha_z) * z_out;
        end

        function [z_out] = Apply_W_z(this, z_in)
            z_out = this.Apply_W_z_Acute(z_in);
            z_out = (1 / this.alpha_z) * z_out;
        end

    end

end
