classdef Synthetic_Test_MD_Interface_Elliptic_Prior < HDSA_Sabl_MD_Interface_Elliptic_Prior

    properties
        con_lofi
        con_hifi
        opt_lofi
        E_u
        M_u
        E_z
        M_z
        E_d
    end

    methods (Access = public)

        function [u_out] = Apply_E_u_Inverse(this, u_in)
            u_out = linsolve(this.E_u, u_in);
        end

        function [u_out] = Apply_E_u_Inverse_Transpose(this, u_in)
            u_out = linsolve(this.E_u', u_in);
        end

        function [u_out] = Apply_M_u(this, u_in)
            u_out = this.M_u * u_in;
        end

        function [u_out] = Apply_M_u_Inverse(this, u_in)
            u_out = linsolve(this.M_u, u_in);
        end

        function [z_out] = Apply_E_z_Inverse(this, z_in)
            z_out = linsolve(this.E_z, z_in);
        end

        function [z_out] = Apply_E_z_Inverse_Transpose(this, z_in)
            z_out = linsolve(this.E_z', z_in);
        end

        function [z_out] = Apply_M_z(this, z_in)
            z_out = this.M_z * z_in;
        end

        function [u_out] = Apply_E_d(this, u_in)
            u_out = this.E_d * u_in;
        end

        function [u_out] = Apply_E_d_Transpose(this, u_in)
            u_out = this.E_d' * u_in;
        end

        function [u_opt] = Load_Optimal_u(this)
            u_opt = load('Opt_Data.mat', 'u').u;
        end

        function [z_opt] = Load_Optimal_z(this)
            z_opt = load('Opt_Data.mat', 'z').z;
        end

        function [Z] = Load_Z_Data(this)
            Z = load('Opt_Data.mat', 'Z').Z;
        end

        function [D] = Load_d_Data(this)
            D = load('Opt_Data.mat', 'D').D;
        end

    end

    methods

        function this = Synthetic_Test_MD_Interface_Elliptic_Prior(opt, alpha_u, alpha_z, con_lofi, con_hifi, opt_lofi)
            this@HDSA_Sabl_MD_Interface_Elliptic_Prior(opt, alpha_u, alpha_z);
            this.con_lofi = con_lofi;
            this.con_hifi = con_hifi;
            this.opt_lofi = opt_lofi;
            m = this.con_lofi.m;
            n = this.con_lofi.n;

            E_u = diag(2 * ones(m, 1)) + diag(-.5 * ones(m - 1, 1), -1) + diag(-.25 * ones(m - 1, 1), 1);
            M_u = eye(m); % diag(1:m);
            E_z = diag(1 * ones(n, 1)) + diag(-.4 * ones(n - 1, 1), -1) + diag(-.2 * ones(n - 1, 1), 1);
            M_z = diag(1:n);
            E_d = .5 * E_u;

            this.E_u = E_u;
            this.M_u = M_u;
            this.E_z = E_z;
            this.M_z = M_z;
            this.E_d = E_d;

            num_sing_vals = 10;
            oversampling = 0;
            num_subspace_iters = 1;
            u_vec = zeros(this.m, 1);
            this.Compute_Elliptic_GSVD(num_sing_vals, oversampling, num_subspace_iters, u_vec);

        end

    end

end
