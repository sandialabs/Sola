classdef Thermochemical_Elliptic_z_Prior_Interface < MD_Elliptic_z_Prior_Interface

    properties
        M
        E_space
        gl_weights
    end

    methods (Access = public)

        function [z_out] = Apply_E_z_Inverse(this, z_in)
            z_out = linsolve(kron(diag(this.gl_weights), this.E_space), z_in);
        end

        function [z_out] = Apply_E_z_Inverse_Transpose(this, z_in)
            z_out = linsolve(kron(diag(this.gl_weights), this.E_space)', z_in);
        end

        function [z_out] = Apply_M_z(this, z_in)
            z_out = kron(diag(this.gl_weights), this.M) * z_in;
        end

        function [z_out] = Apply_E_z(this, z_in)
            z_out = kron(diag(this.gl_weights), this.E_space) * z_in;
        end

        function [z_out] = Apply_E_z_Transpose(this, z_in)
            z_out = kron(diag(this.gl_weights), this.E_space)' * z_in;
        end

        function [z_out] = Apply_M_z_Inverse(this, z_in)
            z_out = linsolve(kron(diag(this.gl_weights), this.M), z_in);
        end

        function this = Thermochemical_Elliptic_z_Prior_Interface(alpha_z, fe, T, control_time_nodes)
            this@MD_Elliptic_z_Prior_Interface(alpha_z);
            this.M = fe.M;
            this.E_space = (1.e-2) * fe.S + fe.M;
            this.gl_weights = (T / 2) * (1 ./ ((0:(control_time_nodes - 1))' + 1 / 2));
        end

    end

end
