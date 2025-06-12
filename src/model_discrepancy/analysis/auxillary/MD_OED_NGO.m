classdef MD_OED_NGO < handle

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

        function this = MD_OED_NGO(opt_prob_interface, data_interface, u_prior_interface, z_prior_interface, md_hessian_analysis, oed_interface)
            arguments
                opt_prob_interface MD_Opt_Prob_Interface
                data_interface MD_Data_Interface
                u_prior_interface MD_u_Prior_Interface
                z_prior_interface MD_z_Prior_Interface
                md_hessian_analysis MD_Hessian_Analysis
                oed_interface MD_OED_Interface
            end
            this.opt_prob_interface = opt_prob_interface;
            this.data_interface = data_interface;
            this.u_prior_interface = u_prior_interface;
            this.z_prior_interface = z_prior_interface;
            this.md_hessian_analysis = md_hessian_analysis;
            this.oed_interface = oed_interface;
            this.verbosity = false;
        end

        function this = Offline_Computation(this)

            this.offline_data = struct;

            this.offline_data.V = this.md_hessian_analysis.evecs;
            this.offline_data.Rho = this.md_hessian_analysis.evals;
            this.offline_data.r = length(this.offline_data.Rho);
            this.offline_data.m = length(this.data_interface.u_opt);

            tmp = this.u_prior_interface.Apply_W_u_Inverse(eye(this.offline_data.m));
            Mu_Wu_inv = this.u_prior_interface.Apply_M_u(tmp);

            V_acute = zeros(length(this.opt_prob_interface.u_current), this.offline_data.r);
            Mu_Wu_inv_V_acute = zeros(length(this.opt_prob_interface.u_current), this.offline_data.r);
            Mz_V = this.z_prior_interface.Apply_M_z(this.offline_data.V);
            Vt_Mz_Wz_inv_Mz_V = zeros(this.offline_data.r, this.offline_data.r);
            Mz_Wz_inv_Mz_V = zeros(this.offline_data.m, this.offline_data.r);
            Wz_inv_Mz_V = zeros(this.offline_data.m, this.offline_data.r);

            for k = 1:this.offline_data.r
                tmp = this.opt_prob_interface.Apply_Solution_Operator_z_Jacobian(this.offline_data.V(:, k), this.data_interface.z_opt);
                V_acute(:, k) = this.opt_prob_interface.Apply_Misfit_Hessian(tmp, this.data_interface.u_opt, this.data_interface.z_opt);

                tmp = this.u_prior_interface.Apply_W_u_Inverse(V_acute(:, k));
                Mu_Wu_inv_V_acute(:, k) = this.u_prior_interface.Apply_M_u(tmp);

                tmp = this.z_prior_interface.Apply_W_z_Inverse(Mz_V(:, k));
                Wz_inv_Mz_V(:, k) = tmp;
                Mz_Wz_inv_Mz_V(:, k) = this.z_prior_interface.Apply_M_z(tmp);
                Vt_Mz_Wz_inv_Mz_V(:, k) = Mz_V' * tmp;
            end

            this.offline_data.V_accute = V_acute;
            this.offline_data.Mu_Wu_inv_V_acute = Mu_Wu_inv_V_acute;
            this.offline_data.Mu_Wu_inv = Mu_Wu_inv;
            this.offline_data.Vt_Mz_Wz_inv_Mz_V = Vt_Mz_Wz_inv_Mz_V;
            this.offline_data.Mz_Wz_inv_Mz_V = Mz_Wz_inv_Mz_V;
            this.offline_data.Wz_inv_Mz_V = Wz_inv_Mz_V;

            Ju = this.opt_prob_interface.Misfit_Gradient(this.data_interface.u_opt, this.data_interface.z_opt);
            tmp = this.u_prior_interface.Apply_W_u_Inverse(Ju);
            this.offline_data.Ju = Ju;
            this.offline_data.Mu_Wu_inv_Ju = this.u_prior_interface.Apply_M_u(tmp);

            this.offline_data.Vt_Design_Cov_Inv_V = this.offline_data.V' * this.oed_interface.Apply_Design_Cov_Inverse(this.offline_data.V);

        end

        function [beta, Z] = Generate_Optimal_Design(this, beta_0, alpha_d, reg_coeff)
            if this.verbosity
                options = optimoptions('fminunc', 'Algorithm', 'quasi-newton', 'Display', 'iter', 'MaxIterations', 5000, 'SpecifyObjectiveGradient', true);
            else
                options = optimoptions('fminunc', 'Algorithm', 'quasi-newton', 'Display', 'None', 'MaxIterations', 5000, 'SpecifyObjectiveGradient', true);
            end
            beta = fminunc(@(beta) this.Evaluate_OED_Objective(beta, alpha_d, reg_coeff), beta_0, options);
            Z = this.data_interface.z_opt + this.offline_data.V * reshape(beta, size(this.offline_data.V, 2), []);
            Z = [this.data_interface.z_opt, Z];
        end

        function [val, grad] = Evaluate_OED_Objective(this, beta, alpha_d, reg_coeff)
            [val1, grad1] = this.Evaluate_Posterior_Cov_Trace(beta, alpha_d);
            [val2, grad2] = this.Evaluate_Regularization(beta);
            val = real(-val1 + reg_coeff * val2);
            grad = real(-grad1 + reg_coeff * grad2);
        end

        function [val, grad] = Evaluate_Posterior_Cov_Trace(this, beta, alpha_d)
            N = length(beta) / this.offline_data.r + 1;

            [g, mu, Mg, g_jac, mu_jac, Mg_jac] = this.G_eigs(beta);

            Ws_Mu_Wu_inv = cell(N, 1); % new
            y_P_y = zeros(N, 1); % new
            s = zeros(N, 1); % new
            p = zeros(N, 1); % new

            for i = 1:N
                Ws_Mu_Wu_inv{i} = (1 / alpha_d) * this.u_prior_interface.Apply_W_u_Plus_scalar_M_u_Inverse(this.offline_data.Mu_Wu_inv, mu(i) / alpha_d); % new
                tmp = this.offline_data.Vt_Mz_Wz_inv_Mz_V * Mg(:, i);
                y_P_y(i) = tmp' * tmp; % new
                s(i) = sum(g(:, i)) - (this.offline_data.Mz_Wz_inv_Mz_V * Mg(:, i))' * this.data_interface.z_opt;
                p(i) = s(i)^2 + y_P_y(i);
            end

            val = 0;
            grad = 0 * beta;
            % Implement gradient via product rules, may be able to
            % precompute additional vectors needed in the analysis
            for i = 1:N
                val = val + p(i) * trace(Ws_Mu_Wu_inv{i});
                grad_si = sum(g_jac{i}, 1)' - Mg_jac{i}' * (this.offline_data.Mz_Wz_inv_Mz_V' * this.data_interface.z_opt);
                grad_yPyi = Mg_jac{i}' * (2 * this.offline_data.Vt_Mz_Wz_inv_Mz_V * (this.offline_data.Vt_Mz_Wz_inv_Mz_V * Mg(:, i)));
                grad_pi = 2 * s(i) * grad_si + grad_yPyi;
                grad = grad + grad_pi * trace(Ws_Mu_Wu_inv{i});

                tmp1 = this.u_prior_interface.Apply_M_u(Ws_Mu_Wu_inv{i});
                tmp2 = -(1 / alpha_d) * this.u_prior_interface.Apply_W_u_Plus_scalar_M_u_Inverse(tmp1, mu(i) / alpha_d);
                grad = grad + p(i) * trace(tmp2) * mu_jac{i};
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
            G = ones(N, N) + M' * this.offline_data.Vt_Mz_Wz_inv_Mz_V * M;
            [g, mu] = eig(G, 'vector');
            Mg = M * g;

            g_jac = cell(N, 1);
            mu_jac = cell(N, 1);
            Mg_jac = cell(N, 1);
            for i = 1:N
                vec2 = this.offline_data.Vt_Mz_Wz_inv_Mz_V * M * g(:, i);
                mu_jac{i} = 2 * kron(g(2:end, i), vec2);

                mat = g * pinv(mu(i) * eye(N) - diag(mu)) * g';
                mat2 = mat * M' * this.offline_data.Vt_Mz_Wz_inv_Mz_V;
                g_jac{i} = kron(mat(:, 2:end), vec2') + kron(g(2:end, i)', mat2);

                Mg_jac{i} = kron(g(2:end, i)', eye(this.offline_data.r)) + M * g_jac{i};
            end

        end

    end

end
