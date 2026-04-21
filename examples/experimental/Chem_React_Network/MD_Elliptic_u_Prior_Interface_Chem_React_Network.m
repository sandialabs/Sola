%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%      Sola - Sandbox for Outer Loop Analysis         %%%%%%%%%
%%%%%%%%% Questions? Contact Joseph Hart (joshart@sandia.gov) %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef MD_Elliptic_u_Prior_Interface_Chem_React_Network < MD_Elliptic_u_Prior_Interface

    properties
        sola_opt
        E_u
        M
    end

    methods

        function [u_out] = Apply_E_u_Inverse(this, u_in)
            u_out = linsolve(this.E_u, u_in);
        end

        function [u_out] = Apply_E_u_Inverse_Transpose(this, u_in)
            u_out = linsolve(this.E_u', u_in);
        end

        function [u_out] = Apply_M_u(this, u_in)
            u_out = kron(eye(this.sola_opt.con.n_y), this.M) * u_in;
        end

        function this = MD_Elliptic_u_Prior_Interface_Chem_React_Network(alpha_u, sola_opt)
            this@MD_Elliptic_u_Prior_Interface(alpha_u);

            this.sola_opt = sola_opt;

            n_t = sola_opt.con.n_t;
            h = sola_opt.con.t_mesh(2) - sola_opt.con.t_mesh(1);
            M = diag(4 * ones(1, n_t)) + diag(ones(1, n_t - 1), 1) + diag(ones(1, n_t - 1), -1);
            M(1, 1) = .5 * M(1, 1);
            M(end, end) = .5 * M(end, end);
            M = (1 / 6) * h * M;
            this.M = M;

            S = diag(2 * ones(1, n_t)) + (-1) * diag(ones(1, n_t - 1), 1) + (-1) * diag(ones(1, n_t - 1), -1);
            S(1, 1) = .5 * S(1, 1);
            S(end, end) = .5 * S(end, end);
            S = (1 / h) * S;

            I = eye(9);
            this.E_u = (5.e-2) * kron(S, I) + kron(this.M, I);

            num_sing_vals = 900;
            oversampling = 0;
            num_subspace_iters = 1;
            u_vec = zeros(sola_opt.con.n_y * n_t, 1);
            this.Compute_E_u_Inverse_GSVD(num_sing_vals, oversampling, num_subspace_iters, u_vec);
        end

    end

end
