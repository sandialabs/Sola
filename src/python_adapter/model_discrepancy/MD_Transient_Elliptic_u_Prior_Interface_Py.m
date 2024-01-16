classdef MD_Transient_Elliptic_u_Prior_Interface_Py < MD_Transient_Elliptic_u_Prior_Interface

    properties
        transient_elliptic_u_prior_interface_py
    end

    methods

        function [u_out] = Apply_E_u_Inverse(this, u_in)
            u_out = this.transient_elliptic_u_prior_interface_py.Apply_E_u_Inverse(u_in);
            u_out = double(u_out);
            if size(u_out, 1) == 1
                u_out = u_out';
            end
        end

        function [u_out] = Apply_E_u_Inverse_Transpose(this, u_in)
            u_out = this.transient_elliptic_u_prior_interface_py.Apply_E_u_Inverse_Transpose(u_in);
            u_out = double(u_out);
            if size(u_out, 1) == 1
                u_out = u_out';
            end
        end

        function [u_out] = Apply_Spatial_M_u(this, u_in)
            u_out = this.transient_elliptic_u_prior_interface_py.Apply_Spatial_M_u(u_in);
            u_out = double(u_out);
            if size(u_out, 1) == 1
                u_out = u_out';
            end
        end

        function this = MD_Transient_Elliptic_u_Prior_Interface_Py(transient_elliptic_u_prior_interface_py, alpha_u, transient_prior_cov)
            this@MD_Transient_Elliptic_u_Prior_Interface(alpha_u, transient_prior_cov);
            this.transient_elliptic_u_prior_interface_py = transient_elliptic_u_prior_interface_py;
        end

    end

end
