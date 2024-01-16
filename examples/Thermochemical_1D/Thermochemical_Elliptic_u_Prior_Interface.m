classdef Thermochemical_Elliptic_u_Prior_Interface < MD_Elliptic_u_Prior_Interface

    properties
        M
        E_T
        E_u1
        E_v1
        E_v2
        m
    end

    methods (Access = public)

        function [u_out] = Apply_E_u_Inverse(this, u_in)
            u_out = 0.0 * u_in;
            u_out(1:this.m, :) = linsolve(this.E_T, u_in(1:this.m, :));
            u_out((this.m + 1):(2 * this.m), :) = linsolve(this.E_u1, u_in((this.m + 1):(2 * this.m), :));
            u_out((2 * this.m + 1):(3 * this.m), :) = linsolve(this.E_v1, u_in((2 * this.m + 1):(3 * this.m), :));
            u_out((3 * this.m + 1):(4 * this.m), :) = linsolve(this.E_v2, u_in((3 * this.m + 1):(4 * this.m), :));
        end

        function [u_out] = Apply_E_u_Inverse_Transpose(this, u_in)
            u_out = 0.0 * u_in;
            u_out(1:this.m, :) = linsolve(this.E_T, u_in(1:this.m, :));
            u_out((this.m + 1):(2 * this.m), :) = linsolve(this.E_u1, u_in((this.m + 1):(2 * this.m), :));
            u_out((2 * this.m + 1):(3 * this.m), :) = linsolve(this.E_v1, u_in((2 * this.m + 1):(3 * this.m), :));
            u_out((3 * this.m + 1):(4 * this.m), :) = linsolve(this.E_v2, u_in((3 * this.m + 1):(4 * this.m), :));
        end

        function [u_out] = Apply_M_u(this, u_in)
            u_out = kron(eye(4), this.M) * u_in;
        end

        function this = Thermochemical_Elliptic_u_Prior_Interface(alpha_u, transient_prior_cov, fe)
            this@MD_Elliptic_u_Prior_Interface(alpha_u, transient_prior_cov);
            this.M = fe.M;
            this.E_T = (5.e-2) * fe.S + fe.M;
            this.E_u1 = (5.e-2) * fe.S + fe.M;
            this.E_v1 = (5.e-2) * fe.S + fe.M;
            this.E_v2 = (5.e-2) * fe.S + fe.M;
            this.m = fe.m;

            num_sing_vals = 200;
            oversampling = 0;
            num_subspace_iters = 1;
            u_vec = zeros(4 * fe.m, 1);
            this.Compute_E_u_Inverse_GSVD(num_sing_vals, oversampling, num_subspace_iters, u_vec);

        end

    end

end
