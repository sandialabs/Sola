classdef MD_Elliptic_u_Prior_Interface_Py < MD_Elliptic_u_Prior_Interface

    properties
        elliptic_u_prior_interface_py
    end

    methods

        function [u_out] = Apply_E_u_Inverse(this, u_in)
            u_out = this.elliptic_u_prior_interface_py.Apply_E_u_Inverse(u_in);
            u_out = double(u_out);
            if size(u_out, 1) == 1
                u_out = u_out';
            end
        end

        function [u_out] = Apply_E_u_Inverse_Transpose(this, u_in)
            u_out = this.elliptic_u_prior_interface_py.Apply_E_u_Inverse_Transpose(u_in);
            u_out = double(u_out);
            if size(u_out, 1) == 1
                u_out = u_out';
            end
        end

        function [u_out] = Apply_M_u(this, u_in)
            u_out = this.elliptic_u_prior_interface_py.Apply_M_u(u_in);
            u_out = double(u_out);
            if size(u_out, 1) == 1
                u_out = u_out';
            end
        end

        function this = MD_Elliptic_u_Prior_Interface_Py(elliptic_u_prior_interface_py, alpha_u)
            this@MD_Elliptic_u_Prior_Interface(alpha_u);
            this.elliptic_u_prior_interface_py = elliptic_u_prior_interface_py;
        end

    end

end
