classdef MD_Scaled_z_Prior_Interface < MD_z_Prior_Interface

    properties
        alpha_z
    end

    methods (Abstract, Access = public)

        %% Pure virtual functions
        
        [z_out] = Apply_W_z_Acute_Inverse(this, z_in)


    end

    methods

        function [] = Set_alpha_z(this,alpha_z_new)
            this.alpha_z = alpha_z_new;
        end

        function [z_out] = Apply_W_z_Inverse(this, z_in)
            z_out = this.alpha_z * this.Apply_W_z_Acute_Inverse(z_in);
        end

        function this = MD_Scaled_z_Prior_Interface(alpha_z)
            this.alpha_z = alpha_z;
        end

    end

end
