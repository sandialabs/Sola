classdef MD_OED < handle

    properties
        opt_prob_interface
        data_interface
        u_prior_interface
        z_prior_interface
        md_hessian_analysis

        offline_data
        verbosity

        covar_coeff
    end

    methods

        function this = MD_OED(opt_prob_interface, data_interface, u_prior_interface, z_prior_interface, md_hessian_analysis)
            arguments
                opt_prob_interface MD_Opt_Prob_Interface
                data_interface MD_Data_Interface
                u_prior_interface MD_u_Prior_Interface
                z_prior_interface MD_z_Prior_Interface
                md_hessian_analysis MD_Hessian_Analysis
            end
            this.opt_prob_interface = opt_prob_interface;
            this.data_interface = data_interface;
            this.u_prior_interface = u_prior_interface;
            this.z_prior_interface = z_prior_interface;
            this.md_hessian_analysis = md_hessian_analysis;
            this.verbosity = false;
            this.covar_coeff = 1;
        end

        function Set_Covariance_Coefficient(this, covar_coeff)
            this.covar_coeff = covar_coeff;
        end

        function this = Offline_Computation(this)

            this.offline_data = struct;
            this.offline_data.V = this.md_hessian_analysis.evecs;
            this.offline_data.r = length(this.md_hessian_analysis.evals);

            Mz_V = this.z_prior_interface.Apply_M_z(this.offline_data.V);
            Wz_inv_Mz_V = this.z_prior_interface.Apply_W_z_Inverse(Mz_V);
            this.offline_data.Mz_Wz_inv_Mz_V = this.z_prior_interface.Apply_M_z(Wz_inv_Mz_V);
            this.offline_data.Vt_Mz_Wz_inv_Mz_V = Mz_V' * Wz_inv_Mz_V;

            try
                sing_vecs = this.u_prior_interface.sing_vecs_output;
                this.offline_data.xjTxj = sum(sing_vecs .* this.u_prior_interface.Apply_M_u(sing_vecs), 1)';
                this.offline_data.lambda_js = 1 ./ (this.u_prior_interface.alpha_u * this.u_prior_interface.sing_vals.^2);
            catch
                sing_vecs = kron(this.u_prior_interface.spatial_prior_cov.sing_vecs_output, this.u_prior_interface.transient_prior_cov.evecs);
                this.offline_data.xjTxj = sum(sing_vecs .* this.u_prior_interface.Apply_M_u(sing_vecs), 1)';
                this.offline_data.lambda_js = 1 ./ (this.u_prior_interface.alpha_u * kron(this.u_prior_interface.spatial_prior_cov.sing_vals.^2, this.u_prior_interface.transient_prior_cov.evals));
            end

        end

        function [beta_new, Z_new] = Generate_Seq_Optimal_Design(this, beta_0, alpha_d, betas, beta_bar, constr_radius)
            if this.verbosity
                options = optimoptions('fmincon', 'Display', 'iter', 'MaxIterations', 5000, 'SpecifyObjectiveGradient', true, 'SpecifyConstraintGradient', true);
            else
                options = optimoptions('fmincon', 'Display', 'None', 'MaxIterations', 5000, 'SpecifyObjectiveGradient', true, 'SpecifyConstraintGradient', true);
            end
            p = length(beta_0) / size(this.offline_data.V, 2);
            fun = @(beta_new) this.Evaluate_OED_Objective_Seq([betas; beta_new], alpha_d, beta_bar, p);
            tmp = @(b) this.offline_data.V' * this.z_prior_interface.Apply_M_z(this.offline_data.V * (b - beta_bar));
            nonlcon = @(b) deal((b - beta_bar)' * tmp(b) - constr_radius, [], 2 * tmp(b), []);
            beta_new = fmincon(fun, beta_0, [], [], [], [], [], [], nonlcon, options);
            Z_new = this.data_interface.z_opt + this.offline_data.V * reshape(beta_new, size(this.offline_data.V, 2), []);
        end

        function [val, grad] = Evaluate_OED_Objective_Seq(this, beta, alpha_d, beta_bar, p)
            [val_full, grad_full] = this.Evaluate_Posterior_Cov_Trace(beta, alpha_d, beta_bar);
            grad = real(-grad_full(end - p * this.offline_data.r + 1:end));
            val = real(-val_full);
        end

        function [val, grad] = Evaluate_Posterior_Cov_Trace(this, beta, alpha_d, beta_bar)
            N = length(beta) / this.offline_data.r + 1;

            [g, mu, Mg, g_jac, mu_jac, Mg_jac] = this.G_eigs(beta);
            tr_Ws_Mu_Wu_inv = zeros(N, 1);
            y_P_y = zeros(N, 1);
            s = zeros(N, 1);
            p = zeros(N, 1);

            for i = 1:N
                tr_Ws_Mu_Wu_inv(i) = sum(this.offline_data.xjTxj ./ (this.offline_data.lambda_js .* (mu(i) +  alpha_d * this.offline_data.lambda_js)));

                tmp = this.offline_data.Mz_Wz_inv_Mz_V * Mg(:, i);
                y_P_y(i) = this.covar_coeff * (tmp' * this.z_prior_interface.Apply_W_z_Inverse(tmp));
                s(i) = sum(g(:, i)) + (beta_bar' * this.offline_data.Vt_Mz_Wz_inv_Mz_V) * Mg(:, i); % Technically, s(i) + z_bar^T (MW^{-1}M) yi
                p(i) = s(i)^2 + y_P_y(i);
            end

            val = 0;
            grad = 0 * beta;
            % Implement gradient via product rules, may be able to
            % precompute additional vectors needed in the analysis
            for i = 1:N
                val = val + p(i) * tr_Ws_Mu_Wu_inv(i);
                grad_si = sum(g_jac{i}, 1)' + Mg_jac{i}' * (this.offline_data.Vt_Mz_Wz_inv_Mz_V * beta_bar);
                tmp = this.covar_coeff * this.offline_data.Mz_Wz_inv_Mz_V' * this.z_prior_interface.Apply_W_z_Inverse(this.offline_data.Mz_Wz_inv_Mz_V * Mg(:, i));
                grad_yPyi = Mg_jac{i}' * (2 * tmp);
                grad_pi = 2 * s(i) * grad_si + grad_yPyi;
                grad = grad + grad_pi * tr_Ws_Mu_Wu_inv(i);

                tmp = -sum(this.offline_data.xjTxj ./ (this.offline_data.lambda_js .* (mu(i) + alpha_d * this.offline_data.lambda_js).^2));
                grad = grad + p(i) * trace(tmp) * mu_jac{i};
            end

        end

        function [g, mu, Mg, g_jac, mu_jac, Mg_jac] = G_eigs(this, beta)
            N = length(beta) / this.offline_data.r + 1;
            M = zeros(this.offline_data.r, N);
            M(:, 2:end) = reshape(beta, this.offline_data.r, N - 1);
            G = ones(N, N) + M' * this.offline_data.Vt_Mz_Wz_inv_Mz_V * M;
            G = (G + G') / 2;
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
