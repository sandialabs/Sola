classdef MD_OED_DeltaCov < handle

    properties
        opt_prob_interface
        data_interface
        u_prior_interface
        z_prior_interface
        md_hessian_analysis
        oed_interface

        offline_data
        verbosity

        covar_coeff
    end

    methods

        function this = MD_OED_DeltaCov(opt_prob_interface, data_interface, u_prior_interface, z_prior_interface, md_hessian_analysis, oed_interface)
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
            this.covar_coeff = 1;
        end

        function Set_Covariance_Coefficient(this, covar_coeff)
            this.covar_coeff = covar_coeff;
        end

        function this = Offline_Computation(this)

            this.offline_data = struct;

            this.offline_data.V = this.md_hessian_analysis.evecs;
            this.offline_data.r = length(this.md_hessian_analysis.evals);
            this.offline_data.m = length(this.data_interface.u_opt);

            Mz_V = this.z_prior_interface.Apply_M_z(this.offline_data.V); % Cost: r matvecs with Mz
            Wz_inv_Mz_V = this.z_prior_interface.Apply_W_z_Inverse(Mz_V); % Cost: r matvecs with Wz_inv

            this.offline_data.Mz_Wz_inv_Mz_V = this.z_prior_interface.Apply_M_z(Wz_inv_Mz_V);
            this.offline_data.Vt_Mz_Wz_inv_Mz_V = Mz_V' * Wz_inv_Mz_V;
            this.offline_data.Vt_Design_Cov_Inv_V = this.offline_data.V' * this.oed_interface.Apply_Design_Cov_Inverse(this.offline_data.V);

            try
                sing_vecs = this.u_prior_interface.sing_vecs_output;
                this.offline_data.xjTxj = sum(sing_vecs .* this.u_prior_interface.Apply_M_u(sing_vecs), 1)';
                this.offline_data.lambda_js = 1 ./ (this.u_prior_interface.alpha_u * this.u_prior_interface.sing_vals.^2);
            catch
                sing_vecs = kron(this.u_prior_interface.spatial_prior_cov.sing_vecs_output, this.u_prior_interface.transient_prior_cov.evecs);
                this.offline_data.xjTxj = sum(sing_vecs .* this.u_prior_interface.Apply_M_u(sing_vecs), 1)';
                this.offline_data.lambda_js = 1 ./ (this.u_prior_interface.alpha_u * kron(this.u_prior_interface.spatial_prior_cov.sing_vals.^2, this.u_prior_interface.transient_prior_cov.evals));
            end

            % W_u_inv = this.u_prior_interface.Apply_M_u(eye(this.offline_data.m));
            % this.offline_data.Mu_Wu_inv = this.u_prior_interface.Apply_M_u(W_u_inv);
        end

        function [beta, Z] = Generate_Optimal_Design(this, beta_0, alpha_d, reg_coeff)
            if this.verbosity
                options = optimoptions('fminunc', 'Algorithm', 'quasi-newton', 'Display', 'iter', 'MaxIterations', 5000, 'SpecifyObjectiveGradient', true);
            else
                options = optimoptions('fminunc', 'Algorithm', 'quasi-newton', 'Display', 'None', 'MaxIterations', 5000, 'SpecifyObjectiveGradient', true);
            end
            beta_bar = zeros(size(this.offline_data.V, 2), 1);
            beta = fminunc(@(beta) this.Evaluate_OED_Objective(beta, alpha_d, reg_coeff, beta_bar), beta_0, options);
            Z = this.data_interface.z_opt + this.offline_data.V * reshape(beta, size(this.offline_data.V, 2), []);
            Z = [this.data_interface.z_opt, Z];
        end

        % NOTE: This does not incorporate sequential prior updates.
        function [beta_new, Z_new] = Generate_Seq_Optimal_Design(this, beta_0, alpha_d, reg_coeff, betas, beta_bar)
            if this.verbosity
                options = optimoptions('fminunc', 'Algorithm', 'quasi-newton', 'Display', 'iter', 'MaxIterations', 5000, 'SpecifyObjectiveGradient', true);
            else
                options = optimoptions('fminunc', 'Algorithm', 'quasi-newton', 'Display', 'None', 'MaxIterations', 5000, 'SpecifyObjectiveGradient', true);
            end
            p = length(beta_0) / size(this.offline_data.V, 2);
            beta_new = fminunc(@(beta_new) this.Evaluate_OED_Objective_Seq([betas; beta_new], alpha_d, reg_coeff, beta_bar, p), beta_0, options);
            Z_new = this.data_interface.z_opt + this.offline_data.V * reshape(beta_new, size(this.offline_data.V, 2), []);
        end

        function [beta_new, Z_new] = Generate_Seq_Optimal_Design_Con_v1(this, beta_0, alpha_d, betas, beta_bar, nonlcon)
            if this.verbosity
                options = optimoptions('fmincon', 'Display', 'iter', 'MaxIterations', 5000, 'SpecifyObjectiveGradient', true, 'SpecifyConstraintGradient', true);
            else
                options = optimoptions('fmincon', 'Display', 'None', 'MaxIterations', 5000, 'SpecifyObjectiveGradient', true, 'SpecifyConstraintGradient', true);
            end
            p = length(beta_0) / size(this.offline_data.V, 2);
            fun = @(beta_new) this.Evaluate_OED_Objective_Seq([betas; beta_new], alpha_d, 0, beta_bar, p);
            beta_new = fmincon(fun, beta_0, [], [], [], [], [], [], nonlcon, options);
            Z_new = this.data_interface.z_opt + this.offline_data.V * reshape(beta_new, size(this.offline_data.V, 2), []);
        end

        function [beta_new, Z_new] = Generate_Seq_Optimal_Design_Con(this, beta_0, alpha_d, reg_coeff, betas, beta_bar, bd_region)
            if this.verbosity
                options = optimoptions('fmincon', 'Display', 'iter', 'MaxIterations', 5000, 'SpecifyObjectiveGradient', true);
            else
                options = optimoptions('fmincon', 'Display', 'None', 'MaxIterations', 5000, 'SpecifyObjectiveGradient', true);
            end
            p = length(beta_0) / size(this.offline_data.V, 2);
            lb = repmat(beta_bar - bd_region, p, 1);
            ub = repmat(beta_bar + bd_region, p, 1);

            fun = @(beta_new) this.Evaluate_OED_Objective_Seq([betas; beta_new], alpha_d, reg_coeff, beta_bar, p);
            beta_new = fmincon(fun, beta_0, [], [], [], [], lb, ub, [], options);
            Z_new = this.data_interface.z_opt + this.offline_data.V * reshape(beta_new, size(this.offline_data.V, 2), []);
        end

        function [beta, Z, post_var, reg_val] = L_Curve_Analysis(this, beta_0, alpha_d, reg_coeffs, betas, beta_bar)
            m = length(reg_coeffs);
            post_var = zeros(m, 1);
            reg_val = zeros(m, 1);
            beta = zeros(length(beta_0), m);
            Z = cell(m, 1);
            for k = 1:m
                [beta(:, k), Z{k}] = this.Generate_Seq_Optimal_Design(beta_0, alpha_d, reg_coeffs(k), betas, beta_bar);
                post_var(k) = -this.Evaluate_Posterior_Cov_Trace(beta(:, k), alpha_d, beta_bar);
                reg_val(k) = this.Evaluate_Regularization(beta(:, k), beta_bar);
            end
        end

        % NOTE: This does not incorporate sequential prior updates.
        function [val, grad] = Evaluate_OED_Objective_Seq(this, beta, alpha_d, reg_coeff, beta_bar, p)
            [val, grad_full] = this.Evaluate_OED_Objective(beta, alpha_d, reg_coeff, beta_bar);
            grad = grad_full(end - p * this.offline_data.r + 1:end);
        end

        function [val, grad] = Evaluate_OED_Objective(this, beta, alpha_d, reg_coeff, beta_bar)
            [val1, grad1] = this.Evaluate_Posterior_Cov_Trace(beta, alpha_d, beta_bar);
            [val2, grad2] = this.Evaluate_Regularization(beta, beta_bar);
            val = real(-val1 + reg_coeff * val2);
            grad = real(-grad1 + reg_coeff * grad2);
        end

        function [val, grad] = Evaluate_Posterior_Cov_Trace(this, beta, alpha_d, beta_bar)
            N = length(beta) / this.offline_data.r + 1;

            [g, mu, Mg, g_jac, mu_jac, Mg_jac] = this.G_eigs(beta);
            %
            Ws_Mu_Wu_inv = cell(N, 1); % old
            tr_Ws_Mu_Wu_inv = zeros(N, 1); % new
            y_P_y = zeros(N, 1); % new
            s = zeros(N, 1); % new
            p = zeros(N, 1); % new

            for i = 1:N
                % Ws_Mu_Wu_inv{i} = (1 / alpha_d) * this.u_prior_interface.Apply_W_u_Plus_scalar_M_u_Inverse(this.offline_data.Mu_Wu_inv, mu(i) / alpha_d); % old
                tr_Ws_Mu_Wu_inv(i) = sum(this.offline_data.xjTxj ./ (this.offline_data.lambda_js .* (mu(i) +  alpha_d * this.offline_data.lambda_js)));

                tmp = this.offline_data.Mz_Wz_inv_Mz_V * Mg(:, i);
                % y_P_y(i) = this.covar_coeff * (tmp' * tmp);
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
                % grad_yPyi = this.covar_coeff * Mg_jac{i}' * (2 * this.offline_data.Mz_Wz_inv_Mz_V' * (this.offline_data.Mz_Wz_inv_Mz_V * Mg(:, i))); % old
                tmp = this.covar_coeff * this.offline_data.Mz_Wz_inv_Mz_V' * this.z_prior_interface.Apply_W_z_Inverse(this.offline_data.Mz_Wz_inv_Mz_V * Mg(:, i));
                grad_yPyi = Mg_jac{i}' * (2 * tmp);
                grad_pi = 2 * s(i) * grad_si + grad_yPyi;
                grad = grad + grad_pi * tr_Ws_Mu_Wu_inv(i);

                tmp = -sum(this.offline_data.xjTxj ./ (this.offline_data.lambda_js .* (mu(i) + alpha_d * this.offline_data.lambda_js).^2));
                grad = grad + p(i) * trace(tmp) * mu_jac{i};
            end

        end

        function [val, grad] = Evaluate_Regularization(this, beta, beta_bar)
            N = length(beta) / this.offline_data.r + 1;
            M = reshape([zeros(this.offline_data.r, 1); beta], this.offline_data.r, N) - beta_bar;
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
