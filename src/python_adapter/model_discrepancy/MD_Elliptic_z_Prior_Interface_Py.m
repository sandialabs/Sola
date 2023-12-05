classdef MD_Elliptic_z_Prior_Interface_Py < MD_Elliptic_z_Prior_Interface

    properties
        elliptic_z_prior_interface_py
    end

    methods

        function [z_out] = Apply_E_z_Inverse(this, z_in)
            z_out = this.elliptic_z_prior_interface_py.Apply_E_z_Inverse(z_in);
            z_out = double(z_out);
            if size(z_out, 1) == 1
                z_out = z_out';
            end
        end

        function [z_out] = Apply_E_z_Inverse_Transpose(this, z_in)
            z_out = this.elliptic_z_prior_interface_py.Apply_E_z_Inverse_Transpose(z_in);
            z_out = double(z_out);
            if size(z_out, 1) == 1
                z_out = z_out';
            end
        end

        function [z_out] = Apply_M_z(this, z_in)
            z_out = this.elliptic_z_prior_interface_py.Apply_M_z(z_in);
            z_out = double(z_out);
            if size(z_out, 1) == 1
                z_out = z_out';
            end
        end

        function [z_out] = Apply_E_z(this, z_in)
            z_out = this.elliptic_z_prior_interface_py.Apply_E_z(z_in);
            z_out = double(z_out);
            if size(z_out, 1) == 1
                z_out = z_out';
            end
        end

        function [z_out] = Apply_E_z_Transpose(this, z_in)
            z_out = this.elliptic_z_prior_interface_py.Apply_E_z_Transpose(z_in);
            z_out = double(z_out);
            if size(z_out, 1) == 1
                z_out = z_out';
            end
        end

        function [z_out] = Apply_M_z_Inverse(this, z_in)
            z_out = this.elliptic_z_prior_interface_py.Apply_M_z_Inverse(z_in);
            z_out = double(z_out);
            if size(z_out, 1) == 1
                z_out = z_out';
            end
        end

        function this = MD_Elliptic_z_Prior_Interface_Py(elliptic_z_prior_interface_py, alpha_z)
            this@MD_Elliptic_z_Prior_Interface(alpha_z);
            this.elliptic_z_prior_interface_py = elliptic_z_prior_interface_py;
        end

    end

end
