classdef MD_u_Prior_Interface_Py < MD_u_Prior_Interface

    properties
        u_prior_interface_py
    end

    methods

        function [u_out] = Apply_M_u(this, u_in)
            u_out = this.u_prior_interface_py.Apply_M_u(u_in);
            u_out = double(u_out);
            if size(u_out, 1) == 1
                u_out = u_out';
            end
        end

        function [u_out] = Apply_W_u_Plus_scalar_M_u_Inverse(this, u_in, beta)
            u_out = this.u_prior_interface_py.Apply_W_u_Plus_scalar_M_u_Inverse(u_in, beta);
            u_out = double(u_out);
            if size(u_out, 1) == 1
                u_out = u_out';
            end
        end

        function [u_out] = Apply_W_u_Inverse(this, u_in)
            u_out = this.u_prior_interface_py.Apply_W_u_Inverse(u_in);
            u_out = double(u_out);
            if size(u_out, 1) == 1
                u_out = u_out';
            end
        end

        function [u_out] = Apply_W_u_Inverse_Factor(this, u_in)
            u_out = this.u_prior_interface_py.Apply_W_u_Inverse_Factor(u_in);
            u_out = double(u_out);
            if size(u_out, 1) == 1
                u_out = u_out';
            end
        end

        function [u_out] = Apply_W_u_Plus_scalar_M_u_Inverse_Factor(this, u_in, scalar)
            u_out = this.u_prior_interface_py.Apply_W_u_Plus_scalar_M_u_Inverse_Factor(u_in, scalar);
            u_out = double(u_out);
            if size(u_out, 1) == 1
                u_out = u_out';
            end
        end

        function this = MD_u_Prior_Interface_Py(u_prior_interface_py)
            this.u_prior_interface_py = u_prior_interface_py;
        end

    end

end
