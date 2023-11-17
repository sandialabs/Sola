classdef MD_Elliptic_z_Prior_Interface < MD_z_Prior_Interface

    properties
        alpha_z
    end

    methods (Abstract, Access = public)

        %% Pure virtual functions

        [z_out] = Apply_E_z_Inverse(this, z_in)

        [z_out] = Apply_E_z_Inverse_Transpose(this, z_in)

        [z_out] = Apply_M_z(this, z_in)

    end

    methods

        % This function must be implemented to enable Hessian GEVP
        function [z_out] = Apply_E_z(this, z_in)
            z_out = [];
            disp('Apply_E_z must be implemented to use Hessian GEVP');
        end

        % This function must be implemented to enable Hessian GEVP
        function [z_out] = Apply_E_z_Transpose(this, z_in)
            z_out = [];
            disp('Apply_E_z_Transpose must be implemented to use Hessian GEVP');
        end

        % This function must be implemented to enable Hessian GEVP
        function [z_out] = Apply_M_z_Inverse(this, z_in)
            z_out = [];
            disp('Apply_M_z_Inverse must be implemented to use Hessian GEVP');
        end

        function this = MD_Elliptic_z_Prior_Interface(alpha_z)
            this.alpha_z = alpha_z;
        end

        % Apply W_z^{-1}
        function [z_out] = Apply_W_z_Inverse(this, z_in)
            tmp1 = this.Apply_E_z_Inverse_Transpose(z_in);
            tmp2 = this.Apply_M_z(tmp1);
            z_out = this.alpha_z * this.Apply_E_z_Inverse(tmp2);
        end

        % Apply W_z
        function [z_out] = Apply_W_z(this, z_in)
            tmp1 = this.Apply_E_z(z_in);
            tmp2 = this.Apply_M_z_Inverse(tmp1);
            z_out = (1 / this.alpha_z) * this.Apply_E_z_Transpose(tmp2);
        end

        % Factorize W_z^{-1} = F*F^T, function gives z_out = F*z_in
        function [z_out] = Apply_W_z_Inverse_Factor(this, z_in)
            mass_mat_sqrt = M_z_Sqrt(this);
            tmp = mass_mat_sqrt.Matrix_Sqrt_Apply(z_in);
            z_out = sqrt(this.alpha_z) * this.Apply_E_z_Inverse(tmp);
        end

    end

end
