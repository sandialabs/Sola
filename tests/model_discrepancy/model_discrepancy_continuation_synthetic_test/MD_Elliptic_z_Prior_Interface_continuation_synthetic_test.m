classdef MD_Elliptic_z_Prior_Interface_continuation_synthetic_test < MD_Elliptic_z_Prior_Interface

    properties
        E_z
        M_z
    end

    methods (Access = public)

        function [z_out] = Apply_E_z_Inverse(this, z_in)
            z_out = linsolve(this.E_z, z_in);
        end

        function [z_out] = Apply_E_z_Inverse_Transpose(this, z_in)
            z_out = linsolve(this.E_z', z_in);
        end

        function [z_out] = Apply_M_z(this, z_in)
            z_out = this.M_z * z_in;
        end

    end

    methods

        function this = MD_Elliptic_z_Prior_Interface_continuation_synthetic_test(alpha_z, n)
            this@MD_Elliptic_z_Prior_Interface(alpha_z);
            E_z = diag(1 * ones(n, 1)) + diag(-.4 * ones(n - 1, 1), -1) + diag(-.2 * ones(n - 1, 1), 1);
            M_z = diag(1:n);
            this.E_z = E_z;
            this.M_z = M_z;
        end

    end

end
