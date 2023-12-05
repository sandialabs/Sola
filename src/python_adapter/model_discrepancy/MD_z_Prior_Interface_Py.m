classdef MD_z_Prior_Interface_Py < MD_z_Prior_Interface

    properties
        z_prior_interface_py
    end

    methods 

        function [z_out] = Apply_W_z_Inverse(this, z_in)
            z_out = this.z_prior_interface_py.Apply_W_z_Inverse(z_in);
            z_out = double(z_out);
            if size(z_out,1)==1
                z_out = z_out';
            end
        end

        function [z_out] = Apply_W_z_Inverse_Factor(this, z_in)
            z_out = this.z_prior_interface_py.Apply_W_z_Inverse_Factor(z_in);
            z_out = double(z_out);
            if size(z_out,1)==1
                z_out = z_out';
            end
        end

        function [z_out] = Apply_W_z(this, z_in)
            z_out = this.z_prior_interface_py.Apply_W_z(z_in);
            z_out = double(z_out);
            if size(z_out,1)==1
                z_out = z_out';
            end
        end

        function this = MD_z_Prior_Interface_Py(z_prior_interface_py)
            this.z_prior_interface_py = z_prior_interface_py;
        end

    end

end
