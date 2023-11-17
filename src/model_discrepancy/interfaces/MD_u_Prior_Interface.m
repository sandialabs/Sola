classdef MD_u_Prior_Interface < handle

    properties

    end

    methods (Abstract, Access = public)

        %% Pure virtual functions

        [u_out] = Apply_W_d(this, u_in)

        [u_out] = Apply_W_u_Plus_scalar_W_d_Inverse(this, u_in, beta)

        [u_out] = Apply_W_u_Inverse(this, u_in)

    end

    methods

        % Factorize W_u^{-1}=F*F^T, function gives u_out=F*u_in
        % This function must be implemented to enable posterior update sampling
        function [u_out] = Apply_W_u_Inverse_Factor(this, u_in)
            u_out = [];
            disp('Apply_W_u_Inverse_Factor must be implemented to use sampling algorithms');
        end

        % Factorize (W_u+scalar*W_d)^{-1}=F*F^T, function gives u_out=F*u_in
        % This function must be implemented to enable posterior update sampling
        function [u_out] = Apply_W_u_Plus_scalar_W_d_Inverse_Factor(this, u_in, scalar)
            u_out = [];
            disp('Apply_W_u_Plus_scalar_W_d_Inverse_Factor must be implemented to use sampling algorithms');
        end

        function this = MD_u_Prior_Interface()

        end

    end

end
