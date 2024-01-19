classdef MD_OED < handle

    properties
        opt_prob_interface
        data_interface
        u_prior_interface
        z_prior_interface
        md_hessian_analysis
        oed_interface

        offline_data
        verbosity
    end

    methods

        function this = MD_OED(opt_prob_interface, data_interface, u_prior_interface, z_prior_interface, md_hessian_analysis, oed_interface)
            this.opt_prob_interface = opt_prob_interface;
            this.data_interface = data_interface;
            this.u_prior_interface = u_prior_interface;
            this.z_prior_interface = z_prior_interface;
            this.md_hessian_analysis = md_hessian_analysis;
            this.oed_interface = oed_interface;
            this.verbosity = true;
        end

        function this = Offline_Computation(this)

            this.offline_data = struct;

            this.offline_data.V = this.md_hessian_analysis.evecs;
            this.offline_data.Rho = this.md_hessian_analysis.evals;
            this.offline_data.r = length(this.offline_data.Rho);

            V_acute = 0 * this.offline_data.V;
            Mu_Wu_inv_V_acute = 0 * this.offline_data.V;
            Vt_Wz_inv_V = zeros(this.offline_data.r, this.offline_data.r);
            for k = 1:this.offline_data.r
                tmp = this.opt_prob_interface.Apply_Solution_Operator_z_Jacobian(this.offline_data.V(:, k), this.data_interface.z_opt);
                V_acute(:, k) = this.opt_prob_interface.Apply_Misfit_Hessian(tmp, this.data_interface.u_opt, this.data_interface.z_opt);

                tmp = this.u_prior_interface.Apply_W_u_Inverse(V_acute(:, k));
                Mu_Wu_inv_V_acute(:, k) = this.u_prior_interface.Apply_M_u(tmp);

                tmp = this.z_prior_interface.Apply_W_z_Inverse(this.offline_data.V(:, k));
                Vt_Wz_inv_V(:, k) = this.offline_data.V' * tmp;
            end
            this.offline_data.V_accute = V_acute;
            this.offline_data.Mu_Wu_inv_V_acute = Mu_Wu_inv_V_acute;
            this.offline_data.Vt_Wz_inv_V = Vt_Wz_inv_V;

            Ju = this.opt_prob_interface.Misfit_Gradient(this.data_interface.u_opt, this.data_interface.z_opt);
            tmp = this.u_prior_interface.Apply_W_u_Inverse(Ju);
            this.offline_data.Ju = Ju;
            this.offline_data.Mu_Wu_inv_Ju = this.u_prior_interface.Apply_M_u(tmp);

            this.offline_data.Vt_Design_Cov_Inv_V = this.offline_data.V' * this.oed_interface.Apply_Design_Cov_Inverse(this.offline_data.V);

        end

        function [Z] = Generate_Random_Design(this, N)
            Omega = randn(length(this.data_interface.z_opt), N - 1);
            Z = this.data_interface.z_opt + this.oed_interface.Apply_Design_Cov_Factor(Omega);
            Z = [this.data_interface.z_opt, Z];
        end

        function [Z] = Generate_Random_Design_from_Subspace(this, N)
            v = randn(size(this.offline_data.V, 2), N - 1);
            coeff = length(this.data_interface.z_opt) / trace(this.offline_data.Vt_Design_Cov_Inv_V);
            Z_tmp = sqrt(coeff) * this.offline_data.V * v;
            Z = [this.data_interface.z_opt, this.data_interface.z_opt + Z_tmp];
        end

        function [beta, Z, post_var, reg_val] = L_Curve_Analysis(this, beta_0, alpha_d, reg_coeffs)
            m = length(reg_coeffs);
            post_var = zeros(m, 1);
            reg_val = zeros(m, 1);
            beta = zeros(length(beta_0), m);
            Z = cell(m, 1);
            for k = 1:m
                [beta(:, k), Z{k}] = this.Generate_Optimal_Design(beta_0, alpha_d, reg_coeffs(k));
                post_var(k) = -this.Evaluate_Posterior_Cov_Trace(beta(:, k), alpha_d);
                reg_val(k) = this.Evaluate_Regularization(beta(:, k));
            end
        end

        function [beta, Z] = Generate_Optimal_Design(this, beta_0, alpha_d, reg_coeff)
            if this.verbosity
                options = optimoptions('fminunc', 'Algorithm', 'quasi-newton', 'Display', 'iter', 'MaxIterations', 5000, 'SpecifyObjectiveGradient', true);
            else
                options = optimoptions('fminunc', 'Algorithm', 'quasi-newton', 'Display', 'None', 'MaxIterations', 5000, 'SpecifyObjectiveGradient', true);
            end
            beta = fminunc(@(beta)this.Evaluate_OED_Objective(beta, alpha_d, reg_coeff), beta_0, options);
            Z = this.data_interface.z_opt + this.offline_data.V * reshape(beta, size(this.offline_data.V, 2), []);
            Z = [this.data_interface.z_opt, Z];
        end

        function [val, grad] = Evaluate_OED_Objective(this, beta, alpha_d, reg_coeff)
            [val1, grad1] = this.Evaluate_Posterior_Cov_Trace(beta, alpha_d);
            [val2, grad2] = this.Evaluate_Regularization(beta);
            val = -val1 + reg_coeff * val2;
            grad = -grad1 + reg_coeff * grad2;
        end

        function [val, grad] = Evaluate_Posterior_Cov_Trace(this, beta, alpha_d)
            N = length(beta) / this.offline_data.r + 1;

            [g, mu, Mg, g_jac, mu_jac, Mg_jac] = this.G_eigs(beta);

            m = length(this.offline_data.Mu_Wu_inv_Ju);
            Ws_V_acute = cell(N, 1);
            Ws_Mu_Wu_inv_Ju = zeros(m, N);
            Quz_y = zeros(m, N);
            y_Qz_y = zeros(N, 1);
            c = zeros(N, 1);

            for i = 1:N
                Ws_V_acute{i} = (1 / alpha_d) * this.u_prior_interface.Apply_W_u_Plus_scalar_M_u_Inverse(this.offline_data.V_accute, mu(i) / alpha_d);
                Ws_Mu_Wu_inv_Ju(:, i) = (1 / alpha_d) * this.u_prior_interface.Apply_W_u_Plus_scalar_M_u_Inverse(this.offline_data.Mu_Wu_inv_Ju, mu(i) / alpha_d);

                tmp = this.offline_data.Vt_Wz_inv_V * Mg(:, i);
                Quz_y(:, i) = this.offline_data.V_accute * diag(1 ./ this.offline_data.Rho.^2) * tmp;
                y_Qz_y(i) = tmp' * diag(1 ./ this.offline_data.Rho.^2) * tmp;

                c(i) = sum(g(:, i));
            end

            val = 0;
            grad = 0 * beta;
            % Implement gradient via product rules, may be able to
            % precompute additional vectors needed in the analysis
            for i = 1:N
                %
                tmp = diag(this.offline_data.Mu_Wu_inv_V_acute' * Ws_V_acute{i})' * (1 ./ this.offline_data.Rho.^2);
                val = val + c(i)^2 * tmp;

                grad = grad + 2 * c(i) * sum(g_jac{i}, 1)' * tmp;

                tmp1 = this.u_prior_interface.Apply_M_u(Ws_V_acute{i});
                tmp2 = (1 / alpha_d) * this.u_prior_interface.Apply_W_u_Plus_scalar_M_u_Inverse(tmp1, mu(i) / alpha_d);
                tmp3 = diag(this.offline_data.Mu_Wu_inv_V_acute' * tmp2)' * (1 ./ this.offline_data.Rho.^2);
                grad = grad - c(i)^2 * tmp3 * mu_jac{i};

                %
                tmp = (Quz_y(:, i)' * Ws_Mu_Wu_inv_Ju(:, i));
                val = val + 2 * c(i) * tmp;

                grad = grad + 2 * sum(g_jac{i}, 1)' * tmp;

                tmp1 = this.u_prior_interface.Apply_M_u(Ws_Mu_Wu_inv_Ju(:, i));
                tmp2 = (1 / alpha_d) * this.u_prior_interface.Apply_W_u_Plus_scalar_M_u_Inverse(tmp1, mu(i) / alpha_d);
                grad = grad - 2 * c(i) * (Quz_y(:, i)' * tmp2) * mu_jac{i};

                tmp3 = this.offline_data.Vt_Wz_inv_V * Mg_jac{i};
                tmp4 = this.offline_data.V_accute * diag(1 ./ this.offline_data.Rho.^2) * tmp3;
                grad = grad + 2 * c(i) * (tmp4' * Ws_Mu_Wu_inv_Ju(:, i));

                %
                tmp = (this.offline_data.Ju' * Ws_Mu_Wu_inv_Ju(:, i));
                val = val + y_Qz_y(i) * tmp;

                grad = grad - y_Qz_y(i) * (this.offline_data.Ju' * tmp2) * mu_jac{i};

                tmp1 = this.offline_data.Vt_Wz_inv_V * Mg_jac{i};
                tmp2 = this.offline_data.Vt_Wz_inv_V * Mg(:, i);
                grad = grad + 2 * tmp * (tmp1' * diag(1 ./ this.offline_data.Rho.^2) * tmp2);
            end

        end

        function [val, grad] = Evaluate_Regularization(this, beta)
            N = length(beta) / this.offline_data.r + 1;
            M = reshape([zeros(this.offline_data.r, 1); beta], this.offline_data.r, N);
            val = 0;
            grad = 0 * beta;
            for i = 2:N
                tmp = this.offline_data.Vt_Design_Cov_Inv_V * M(:, i);
                val = val + M(:, i)' * tmp;
                ei = zeros(N, 1);
                ei(i) = 1;
                grad_reg = 2 * kron(ei(2:end), tmp);
                grad = grad + grad_reg;
            end
        end

        function [g, mu, Mg, g_jac, mu_jac, Mg_jac] = G_eigs(this, beta)
            N = length(beta) / this.offline_data.r + 1;
            M = zeros(this.offline_data.r, N);
            M(:, 2:end) = reshape(beta, this.offline_data.r, N - 1);
            G = ones(N, N) + M' * this.offline_data.Vt_Wz_inv_V * M;
            [g, mu] = eig(G, 'vector');
            Mg = M * g;

            g_jac = cell(N, 1);
            mu_jac = cell(N, 1);
            Mg_jac = cell(N, 1);
            for i = 1:N
                vec2 = this.offline_data.Vt_Wz_inv_V * M * g(:, i);
                mu_jac{i} = 2 * kron(g(2:end, i), vec2);

                mat = g * pinv(mu(i) * eye(N) - diag(mu)) * g';
                mat2 = mat * M' * this.offline_data.Vt_Wz_inv_V;
                g_jac{i} = kron(mat(:, 2:end), vec2') + kron(g(2:end, i)', mat2);

                Mg_jac{i} = kron(g(2:end, i)', eye(this.offline_data.r)) + M * g_jac{i};
            end

        end

    end

end
