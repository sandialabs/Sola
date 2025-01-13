classdef MD_Scaled_u_Prior_Interface < MD_u_Prior_Interface

    properties
        alpha_u
    end

    methods (Abstract, Access = public)

        %% Pure virtual functions

        [u_out] = Apply_M_u(this, u_in)

        [u_out] = Apply_W_u_Acute_Plus_scalar_M_u_Inverse(this, u_in, scalar)

        [u_out] = Apply_W_u_Acute_Inverse(this, u_in)

    end

    methods

        function [] = Set_alpha_u(this,alpha_u_new)
            this.alpha_u = alpha_u_new;
        end

        function [u_out] = Apply_W_u_Plus_scalar_M_u_Inverse(this, u_in, scalar)
            u_out = this.alpha_u * this.Apply_W_u_Acute_Plus_scalar_M_u_Inverse(u_in, scalar * this.alpha_u);
        end

        function [u_out] = Apply_W_u_Inverse(this, u_in)
            u_out = this.alpha_u * this.Apply_W_u_Acute_Inverse(u_in);
        end

        function this = MD_Scaled_u_Prior_Interface(alpha_u)
            this.alpha_u = alpha_u;
        end

    end

end
