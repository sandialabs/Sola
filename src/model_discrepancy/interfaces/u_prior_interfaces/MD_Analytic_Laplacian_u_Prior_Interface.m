classdef MD_Analytic_Laplacian_u_Prior_Interface < MD_Elliptic_u_Prior_Interface

    properties
        M
        u_hyperparam_interface
        determine_u_hyperparams
        beta_u
    end

    methods (Access = public)

        function [u_out] = Apply_M_u(this, u_in)
            u_out = this.M * u_in;
        end

        function [u_out] = Apply_E_u_Inverse(this, u_in)
            u_out = [];
            disp('Apply_E_u_Inverse should not be called from MD_Analytic_Laplacian_u_Prior_Interface');
        end

        function [u_out] = Apply_E_u_Inverse_Transpose(this, u_in)
            u_out = [];
            disp('Apply_E_u_Inverse_Transpose should not be called from MD_Analytic_Laplacian_u_Prior_Interface');
        end

        function [] = Compute_E_u_Inverse_GSVD(this, num_sing_vals, oversampling, num_subspace_iters, u_vec)
            disp('Compute_E_u_Inverse_GSVD should not be called from MD_Analytic_Laplacian_u_Prior_Interface');
        end

    end

    methods

        function [] = Generate_Spectral_Decomposition(this)

            nodes = this.u_hyperparam_interface.Load_Spatial_Node_Data();
            nodes = nodes{this.u_hyperparam_interface.component_id};
            m = size(nodes, 1);
            d = size(nodes, 2);
            r = this.u_hyperparam_interface.gsvd_num_sing_vals;

            if d == 1

                L = max(nodes(:, 1)) - min(nodes(:, 1));
                x = nodes(:, 1) - min(nodes(:, 1));
                e = 1 + this.beta_u * (pi / L)^2 * (0:(r - 1)).^2;
                this.sing_vals = 1 ./ e';
                this.sing_vecs_output = zeros(m, r);
                v = ones(m, 1);
                this.sing_vecs_output(:, 1) = v / sqrt(v' * this.Apply_M_u(v));
                for k = 2:r
                    v = cos((k - 1) * (pi / L) * x);
                    this.sing_vecs_output(:, k) = v / sqrt(v' * this.Apply_M_u(v));
                end

            elseif d == 2

                rsq = floor(sqrt(r));
                r = zeros(2, 1);
                L = zeros(2, 1);
                for j = 1:2
                    L(j) = max(nodes(:, j)) - min(nodes(:, j));
                    tmp = (L(j) / pi)^2 * (1 / this.beta_u) * (1 / this.u_hyperparam_interface.W_u_inv_spectral_gap - 1);
                    r(j) = round(sqrt(tmp));
                    r(j) = min(r(j), rsq);
                end

                I = (0:(r(1) - 1))';
                evali = kron(I.^2, ones(r(2), 1));
                I = (0:(r(2) - 1))';
                evalj = kron(ones(r(1), 1), I.^2);
                this.sing_vals = 1 ./ (1 + this.beta_u * pi^2 * (evali / L(1)^2 + evalj / L(2)^2));

                x = nodes(:, 1) - min(nodes(:, 1));
                y = nodes(:, 2) - min(nodes(:, 2));
                v1 = cos(pi * x * (0:(r(1) - 1)) / L(1));
                v2 = cos(pi * y * (0:(r(2) - 1)) / L(2));
                this.sing_vecs_output = zeros(m, r(1) * r(2));
                for i = 1:r(1)
                    I = (1 + (i - 1) * r(2)):(i * r(2));
                    this.sing_vecs_output(:, I) = v2 .* (v1(:, i) * ones(1, r(1)));
                end

                [~, I] = sort(this.sing_vals, 'descend');
                this.sing_vals = this.sing_vals(I);
                this.sing_vecs_output = this.sing_vecs_output(:, I);
                normalize = sqrt(this.sing_vecs_output' * this.Apply_M_u(this.sing_vecs_output));
                for l = 1:r(1) * r(2)
                    this.sing_vecs_output(:, l) = this.sing_vecs_output(:, l) / normalize(l, l);
                end

            elseif d == 3

                disp('MD_Analytic_Laplacian_u_Prior_Interface is not supported in 3D');

            end

        end

        function [] = Set_beta_u(this, beta_u_new)
            this.beta_u = beta_u_new;
            this.Generate_Spectral_Decomposition();
        end

        function this = MD_Analytic_Laplacian_u_Prior_Interface(M, data_interface, u_hyperparam_interface)
            this@MD_Elliptic_u_Prior_Interface(u_hyperparam_interface.alpha_u);
            this.M = M;
            this.u_hyperparam_interface = u_hyperparam_interface;
            this.determine_u_hyperparams = MD_Determine_u_Hyperparameters(data_interface, u_hyperparam_interface);
            this.beta_u = 0.0;

            if this.u_hyperparam_interface.beta_u == 0.0
                this.determine_u_hyperparams.Determine_beta_u();
            end
            this.Set_beta_u(this.u_hyperparam_interface.beta_u);

            if this.u_hyperparam_interface.gsvd_num_sing_vals == 0
                this.determine_u_hyperparams.Determine_GSVD_Hyperparameters();
            end

            this.Generate_Spectral_Decomposition();

            if ~this.u_hyperparam_interface.is_transient
                if this.u_hyperparam_interface.alpha_u == 0.0
                    this.determine_u_hyperparams.Determine_alpha_u(this);
                end
                this.Set_alpha_u(this.u_hyperparam_interface.alpha_u);
            end
        end

    end

end
