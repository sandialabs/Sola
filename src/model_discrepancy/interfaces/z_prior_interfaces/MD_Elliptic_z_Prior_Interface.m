classdef MD_Elliptic_z_Prior_Interface < MD_Scaled_z_Prior_Interface

    properties

    end

    methods (Abstract, Access = public)

        %% Pure virtual functions

        [z_out] = Apply_E_z_Inverse(this, z_in)

        [z_out] = Apply_E_z_Inverse_Transpose(this, z_in)

        [z_out] = Apply_M_z(this, z_in)

    end

    methods

        % Apply W_z^{-1}
        function [z_out] = Apply_W_z_Acute_Inverse(this, z_in)
            tmp1 = this.Apply_E_z_Inverse_Transpose(z_in);
            tmp2 = this.Apply_M_z(tmp1);
            z_out = this.Apply_E_z_Inverse(tmp2);
        end

        % This function must be implemented to enable Hessian GEVP
        function [z_out] = Apply_E_z(this, z_in)
            z_out = [];
            disp('MD_Elliptic_z_Prior_Interface::Apply_E_z must be implemented to use Hessian GEVP');
        end

        % This function must be implemented to enable Hessian GEVP
        function [z_out] = Apply_E_z_Transpose(this, z_in)
            z_out = [];
            disp('MD_Elliptic_z_Prior_Interface::Apply_E_z_Transpose must be implemented to use Hessian GEVP');
        end

        % This function must be implemented to enable Hessian GEVP
        function [z_out] = Apply_M_z_Inverse(this, z_in)
            z_out = [];
            disp('MD_Elliptic_z_Prior_Interface::Apply_M_z_Inverse must be implemented to use Hessian GEVP');
        end

        % Apply W_z_Acute matrix
        % This function must be implemented to enable Hessian GEVP
        function [z_out] = Apply_W_z_Acute(this, z_in)
            tmp = this.Apply_E_z(z_in);
            tmp = this.Apply_M_z_Inverse(tmp);
            z_out = this.Apply_E_z_Transpose(tmp);
        end

        function [] = Update_alpha_z(this, alpha_z_new)
            this.alpha_z = alpha_z_new;
        end

        function this = MD_Elliptic_z_Prior_Interface(alpha_z)
            this@MD_Scaled_z_Prior_Interface(alpha_z);
        end

    end

end
