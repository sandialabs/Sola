classdef MD_Elliptic_z_Prior_Interface < MD_Scaled_z_Prior_Interface

    %%%%%%%%%%%%%%%%% Defined covariance as a squared inverse elliptic operator %%%%%%%%%%%%%%%%%

    properties

    end

    %% Pure virtual functions for user implementation
    methods (Abstract, Access = public)

        [z_out] = Apply_M_z(this, z_in)

        [z_out] = Apply_E_z_Inverse(this, z_in)

        [z_out] = Apply_E_z_Inverse_Transpose(this, z_in)

    end

    %% Virtual functions for user implementation
    methods

        function [z_out] = Apply_E_z(this, z_in)
            z_out = [];
            disp('MD_Elliptic_z_Prior_Interface::Apply_E_z must be implemented to use Hessian GEVP');
        end

        function [z_out] = Apply_E_z_Transpose(this, z_in)
            z_out = [];
            disp('MD_Elliptic_z_Prior_Interface::Apply_E_z_Transpose must be implemented to use Hessian GEVP');
        end

        function [z_out] = Sample_with_Covariance_W_z_Acute_Inverse(this, num_samples)
            disp('MD_Elliptic_z_Prior_Interface::Sample_with_Covariance_W_z_Acute_Inverse must be implemented to use sampling algorithms');
        end

        function [z_out] = Apply_M_z_Inverse(this, z_in)
            z_out = [];
            disp('MD_Elliptic_z_Prior_Interface::Apply_M_z_Inverse must be implemented to use Hessian GEVP');
        end

    end

    %% Constructor and helper functions
    methods

        function this = MD_Elliptic_z_Prior_Interface(alpha_z)
            arguments
                alpha_z (1, 1) {mustBeNumeric}
            end
            this@MD_Scaled_z_Prior_Interface(alpha_z);
        end

        function [] = Update_alpha_z(this, alpha_z_new)
            this.alpha_z = alpha_z_new;
        end

    end

    %% Implementation of base class virtual functions
    methods

        function [z_out] = Apply_W_z_Acute_Inverse(this, z_in)
            tmp1 = this.Apply_E_z_Inverse_Transpose(z_in);
            tmp2 = this.Apply_M_z(tmp1);
            z_out = this.Apply_E_z_Inverse(tmp2);
        end

        function [z_out] = Apply_W_z_Acute(this, z_in)
            tmp = this.Apply_E_z(z_in);
            tmp = this.Apply_M_z_Inverse(tmp);
            z_out = this.Apply_E_z_Transpose(tmp);
        end

    end

end
