classdef Mass_Spring_HDSA < HDSA_Sabl_MD_Interface_Elliptic_Prior

    properties
        E_u
        E_z
        E_d
        P_z
        M
        H
        evecs
        evals
    end

    methods

        function this = Mass_Spring_HDSA(con_opt_obj, alpha_u, alpha_z)
            this@HDSA_Sabl_MD_Interface_Elliptic_Prior(con_opt_obj, alpha_u, alpha_z);

            n_t = con_opt_obj.con.n_t;
            h = con_opt_obj.con.t_mesh(2) - con_opt_obj.con.t_mesh(1);
            M = diag(4 * ones(1, n_t)) + diag(ones(1, n_t - 1), 1) + diag(ones(1, n_t - 1), -1);
            M(1, 1) = .5 * M(1, 1);
            M(end, end) = .5 * M(end, end);
            M = (1 / 6) * h * M;
            this.M = M;

            S = diag(2 * ones(1, n_t)) + (-1) * diag(ones(1, n_t - 1), 1) + (-1) * diag(ones(1, n_t - 1), -1);
            S(1, 1) = .5 * S(1, 1);
            S(end, end) = .5 * S(end, end);
            S = (1 / h) * S;

            this.P_z = con_opt_obj.con.P_z;

            I = eye(2);
            I(2, 2) = 4;
            this.E_u = (5.e-2) * kron(S, I) + kron(this.M, I);
            this.E_u(1:2, :) = 0;
            this.E_u(1, 1) = 10;
            this.E_u(2, 2) = 40;

            this.E_z = (1.e-1) * S + this.M;

            this.E_d = (1.e-8) * kron(S, I) + kron(this.M, I);

            num_sing_vals = 100;
            oversampling = 0;
            num_subspace_iters = 1;
            u_vec = zeros(con_opt_obj.this.n_y * n_t, 1);
            this.Compute_Elliptic_GSVD(num_sing_vals, oversampling, num_subspace_iters, u_vec);
        end

        function [u_out] = Apply_E_u_Inverse(this, u_in)
            u_out = linsolve(this.E_u, u_in);
        end

        function [u_out] = Apply_E_u_Inverse_Transpose(this, u_in)
            u_out = linsolve(this.E_u', u_in);
        end

        function [u_out] = Apply_M_u(this, u_in)
            u_out = kron(eye(this.con_opt_obj.con.n_y), this.M) * u_in;
        end

        function [u_out] = Apply_M_u_Inverse(this, u_in)
            u_out = linsolve(kron(eye(this.con_opt_obj.con.n_y), this.M), u_in);
        end

        function [z_out] = Apply_E_z_Inverse(this, z_in)
            z_out = linsolve(this.P_z' * this.E_z * this.P_z, z_in);
        end

        function [z_out] = Apply_E_z_Inverse_Transpose(this, z_in)
            z_out = linsolve(this.P_z' * this.E_z' * this.P_z, z_in);
        end

        function [z_out] = Apply_M_z(this, z_in)
            z_out = this.P_z' * this.M * this.P_z * z_in;
        end

        % This function must be implemented to enable Hessian GEVP
        function [z_out] = Apply_E_z(this, z_in)
            z_out = this.P_z' * this.E_z * this.P_z * z_in;
        end

        % This function must be implemented to enable Hessian GEVP
        function [z_out] = Apply_E_z_Transpose(this, z_in)
            z_out = this.P_z' * this.E_z' * this.P_z * z_in;
        end

        % This function must be implemented to enable Hessian GEVP
        function [z_out] = Apply_M_z_Inverse(this, z_in)
            z_out = linsolve(this.P_z' * this.M * this.P_z, z_in);
        end

        function [u_out] = Apply_E_d(this, u_in)
            u_out = this.E_d * u_in;
        end

        function [u_out] = Apply_E_d_Transpose(this, u_in)
            u_out = this.E_d' * u_in;
        end

        function [u_opt] = Load_Optimal_u(this)
            u_opt = load('Optimization_Results.mat').u_lofi;
        end

        function [z_opt] = Load_Optimal_z(this)
            z_opt = load('Optimization_Results.mat').z_lofi;
        end

        function [Z] = Load_Z_Data(this)
            Z = load('Optimization_Results.mat').Z;
        end

        function [D] = Load_d_Data(this)
            D = load('Optimization_Results.mat').D;
        end

    end

end
