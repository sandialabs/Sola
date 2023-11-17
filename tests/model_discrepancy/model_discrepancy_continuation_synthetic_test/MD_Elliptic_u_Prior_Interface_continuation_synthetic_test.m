classdef MD_Elliptic_u_Prior_Interface_continuation_synthetic_test < MD_Elliptic_u_Prior_Interface

    properties
        E_u
        M_u
        E_d
    end

    methods (Access = public)

        function [u_out] = Apply_E_u_Inverse(this, u_in)
            u_out = linsolve(this.E_u, u_in);
        end

        function [u_out] = Apply_E_u_Inverse_Transpose(this, u_in)
            u_out = linsolve(this.E_u, u_in);
        end

        function [u_out] = Apply_M_u(this, u_in)
            u_out = this.M_u * u_in;
        end

        function [u_out] = Apply_M_u_Inverse(this, u_in)
            u_out = linsolve(this.M_u, u_in);
        end

        function [u_out] = Apply_E_d(this, u_in)
            u_out = this.E_d * u_in;
        end

        function [u_out] = Apply_E_d_Transpose(this, u_in)
            u_out = this.E_d' * u_in;
        end

    end

    methods

        function this = MD_Elliptic_u_Prior_Interface_continuation_synthetic_test(alpha_u, m)
            this@MD_Elliptic_u_Prior_Interface(alpha_u);

            E_u = diag(2 * ones(m, 1)) + diag(-.5 * ones(m - 1, 1), -1) + diag(-.25 * ones(m - 1, 1), 1);
            M_u = eye(m);
            E_d = .5 * E_u;

            this.E_u = E_u;
            this.M_u = M_u;
            this.E_d = E_d;

            num_sing_vals = 10;
            oversampling = 0;
            num_subspace_iters = 1;
            u_vec = zeros(m, 1);
            this.Compute_E_u_Inverse_GSVD(num_sing_vals, oversampling, num_subspace_iters, u_vec);

        end

    end

end
